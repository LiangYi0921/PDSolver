#pragma once

#include "MstModel.h"
#include <boost/spirit/include/qi.hpp>
#include <fstream>
#include <string>
#include <iostream>

namespace MST {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

class Parser {
public:
    static bool parseFile(const std::string& filepath, Model& model) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file: " << filepath << std::endl;
            return false;
        }

        std::string line;
        std::string currentSection;
        size_t lineNum = 0;

        while (std::getline(file, line)) {
            ++lineNum;

            // 去除首尾空格
            size_t start = line.find_first_not_of(" \t\r\n");
            if (start == std::string::npos) continue;
            size_t end = line.find_last_not_of(" \t\r\n");
            std::string trimmed = line.substr(start, end - start + 1);

            // 忽略注释与分隔线
            if (trimmed.empty() || trimmed == "##" || trimmed[0] == '#') {
                continue;
            }

            // 识别段落头
            if (trimmed.front() == '/' && trimmed.back() == '/' && trimmed.length() > 2) {
                currentSection = trimmed.substr(1, trimmed.length() - 2);
                continue;
            }

            auto it = trimmed.begin();
            auto itEnd = trimmed.end();

            if (currentSection == "NODECORD2") {
                Node node;
                bool ok = qi::phrase_parse(
                    it, itEnd,
                    qi::int_ >> qi::double_ >> qi::double_ >> qi::double_ >> qi::int_ >> qi::int_,
                    ascii::blank,
                    node.id, node.x, node.y, node.z, node.layer, node.flag
                );
                if (ok) {
                    node.z = -node.z; // 恢复正常工程右手坐标系：Z轴向上为正
                    model.nodes.push_back(node);
                }
            }
            else if (currentSection == "MAST_ELEMS_6") {
                Element elem;
                bool ok = qi::phrase_parse(
                    it, itEnd,
                    qi::int_ >> qi::int_ >> qi::int_ >> qi::int_ >> qi::int_ >> qi::int_,
                    ascii::blank,
                    elem.id, elem.memberId, elem.n1, elem.n2, elem.secType, elem.secId
                );
                if (ok) model.elements.push_back(elem);
            }
            else if (currentSection == "SECTION PROPERTIES") {
                SectionProperty sp;
                std::vector<double> dvals;
                std::vector<int> ivals;
                bool ok = qi::phrase_parse(
                    it, itEnd,
                    qi::int_ >> qi::int_ >> *qi::double_ >> *qi::int_,
                    ascii::blank,
                    sp.id, sp.secId, dvals, ivals
                );
                if (ok && dvals.size() >= 2) {
                    sp.D = dvals[0];
                    sp.t = dvals[1];
                    if (dvals.size() > 2) {
                        sp.extras.assign(dvals.begin() + 2, dvals.end());
                    }
                    if (ivals.size() >= 1) sp.flag1 = ivals[0];
                    if (ivals.size() >= 2) sp.flag2 = ivals[1];
                    model.sectionProperties.push_back(sp);
                }
            }
            else if (currentSection == "SUPPORT INFORMATION") {
                SupportInfo si;
                bool ok = qi::phrase_parse(
                    it, itEnd,
                    qi::int_ >> qi::int_ >> qi::int_ >> qi::int_ >> qi::int_ >> qi::int_ >> qi::int_,
                    ascii::blank,
                    si.nodeId, si.dx, si.dy, si.dz, si.rx, si.ry, si.flag
                );
                if (ok) model.supportInfos.push_back(si);
            }
            else if (currentSection == "SPRING INFORMATION") {
                std::vector<double> vals;
                qi::phrase_parse(it, itEnd, *qi::double_, ascii::blank, vals);
                model.springInfoParams = vals;
            }
            else if (currentSection == "SPRING SUPPORT") {
                SpringSupport sp;
                bool ok = qi::phrase_parse(
                    it, itEnd,
                    qi::int_ >> qi::double_ >> qi::double_ >> qi::double_,
                    ascii::blank,
                    sp.nodeId, sp.kx, sp.ky, sp.kz
                );
                if (ok) model.springSupports.push_back(sp);
            }
            else if (currentSection == "AREAINFO1") {
                AreaInfo ai;
                std::vector<double> dvals;
                bool ok = qi::phrase_parse(
                    it, itEnd,
                    qi::int_ >> qi::int_ >> qi::int_ >> *qi::double_,
                    ascii::blank,
                    ai.id, ai.secId, ai.flag, dvals
                );
                if (ok && dvals.size() >= 8) {
                    ai.A   = dvals[0];
                    ai.Iy  = dvals[1];
                    ai.Iz  = dvals[2];
                    ai.It  = dvals[3];
                    ai.E   = dvals[4];
                    ai.G   = dvals[5];
                    ai.nu  = dvals[6];
                    ai.rho = dvals[7];
                    model.areaInfos.push_back(ai);
                }
            }
            else if (currentSection == "NODEMASS") {
                NodeMass nm;
                bool ok = qi::phrase_parse(
                    it, itEnd,
                    qi::int_ >> qi::double_,
                    ascii::blank,
                    nm.nodeId, nm.mass
                );
                if (ok) model.nodeMasses.push_back(nm);
            }
            else if (currentSection == "DEADLOAD") {
                NodeLoad nl;
                bool ok = qi::phrase_parse(
                    it, itEnd,
                    qi::int_ >> qi::double_ >> qi::double_ >> qi::double_,
                    ascii::blank,
                    nl.nodeId, nl.fx, nl.fy, nl.fz
                );
                if (ok) {
                    nl.fz = -nl.fz; // 恢复正常方向：重力向下为负
                    model.deadLoads.push_back(nl);
                }
            }
            else if (currentSection == "LIVELOAD") {
                CaseNodeLoad cnl;
                bool ok = qi::phrase_parse(
                    it, itEnd,
                    qi::int_ >> qi::int_ >> qi::double_ >> qi::double_ >> qi::double_,
                    ascii::blank,
                    cnl.caseId, cnl.nodeId, cnl.fx, cnl.fy, cnl.fz
                );
                if (ok) {
                    cnl.fz = -cnl.fz; // 恢复正常方向：向下为负
                    model.liveLoads.push_back(cnl);
                }
            }
            else if (currentSection == "WINDLOAD") {
                CaseNodeLoad cnl;
                bool ok = qi::phrase_parse(
                    it, itEnd,
                    qi::int_ >> qi::int_ >> qi::double_ >> qi::double_ >> qi::double_,
                    ascii::blank,
                    cnl.caseId, cnl.nodeId, cnl.fx, cnl.fy, cnl.fz
                );
                if (ok) {
                    cnl.fz = -cnl.fz; // 恢复正常方向：向上风吸力为正
                    model.windLoads.push_back(cnl);
                }
            }
            else if (currentSection == "END") {
                break;
            }
        }

        model.buildIndex();
        return true;
    }
};

} // namespace MST
