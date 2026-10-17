#pragma once

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace MST {

// 节点数据
struct Node {
    int id = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    int layer = 0; // 1: 上弦层, 2: 下弦层
    int flag = 0;
};

// 单元/杆件数据
struct Element {
    int id = 0;
    int memberId = 0;
    int n1 = 0;
    int n2 = 0;
    int secType = 1; // 倒数第二列: 截面类型 (1 = 圆钢管)
    int secId = 1;   // 最后一列: 截面序号 (1 ~ 8)
};

// 截面几何参数
struct SectionProperty {
    int id = 0;
    int secId = 0;
    double D = 0.0; // 钢管外径 (mm)
    double t = 0.0; // 壁厚 (mm)
    std::vector<double> extras; // 预留参数
    int flag1 = 0;
    int flag2 = 0;
};

// 截面力学与材料特性
struct AreaInfo {
    int id = 0;
    int secId = 0;
    int flag = 0;
    double A = 0.0;     // 截面积 (mm²)
    double Iy = 0.0;    // y轴惯性矩 (mm⁴)
    double Iz = 0.0;    // z轴惯性矩 (mm⁴)
    double It = 0.0;    // 扭转惯性矩 (mm⁴)
    double E = 206000.0;// 弹性模量 (MPa)
    double G = 79000.0; // 剪切模量 (MPa)
    double nu = 0.3;    // 泊松比
    double rho = 7.85e-6; // 密度 (kg/mm³), 对应 7850 kg/m³
};

// 弹性支座 (/SPRING SUPPORT/)
struct SpringSupport {
    int nodeId = 0;
    double kx = 0.0; // kN/mm
    double ky = 0.0; // kN/mm
    double kz = 0.0; // kN/mm
};

// 支座约束信息 (/SUPPORT INFORMATION/)
// dx, dy, dz, rx, ry, flag
// 0 = 自由 (Free), 1 = 刚性固定 (Fixed), 2 = 弹性支承 (Spring)
struct SupportInfo {
    int nodeId = 0;
    int dx = 0;
    int dy = 0;
    int dz = 0;
    int rx = 0;
    int ry = 0;
    int flag = 0;
};

// 节点质量
struct NodeMass {
    int nodeId = 0;
    double mass = 0.0;
};

// 节点单工况荷载 (如恒载)
struct NodeLoad {
    int nodeId = 0;
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
};

// 多工况节点荷载 (如活载、风载)
struct CaseNodeLoad {
    int caseId = 0;
    int nodeId = 0;
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
};

// 构件分类类型
enum class MemberType {
    TopChord,    // 上弦杆
    BottomChord, // 下弦杆
    WebMember    // 腹杆
};

// MSTCAD 模型整体结构
struct Model {
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<SectionProperty> sectionProperties;
    std::vector<AreaInfo> areaInfos;
    std::vector<double> springInfoParams;
    std::vector<SpringSupport> springSupports;
    std::vector<SupportInfo> supportInfos;
    std::vector<NodeMass> nodeMasses;
    std::vector<NodeLoad> deadLoads;
    std::vector<CaseNodeLoad> liveLoads;
    std::vector<CaseNodeLoad> windLoads;

    // 辅助查询缓存
    std::map<int, size_t> nodeIndexMap;
    std::map<int, size_t> sectionIndexMap;
    std::map<int, size_t> areaIndexMap;
    std::map<int, size_t> supportInfoMap;

    void buildIndex() {
        nodeIndexMap.clear();
        for (size_t i = 0; i < nodes.size(); ++i) {
            nodeIndexMap[nodes[i].id] = i;
        }

        sectionIndexMap.clear();
        for (size_t i = 0; i < sectionProperties.size(); ++i) {
            sectionIndexMap[sectionProperties[i].secId] = i;
        }

        areaIndexMap.clear();
        for (size_t i = 0; i < areaInfos.size(); ++i) {
            areaIndexMap[areaInfos[i].secId] = i;
        }

        supportInfoMap.clear();
        for (size_t i = 0; i < supportInfos.size(); ++i) {
            supportInfoMap[supportInfos[i].nodeId] = i;
        }
    }

    const Node* getNode(int id) const {
        auto it = nodeIndexMap.find(id);
        if (it != nodeIndexMap.end()) return &nodes[it->second];
        return nullptr;
    }

    const SectionProperty* getSection(int secId) const {
        auto it = sectionIndexMap.find(secId);
        if (it != sectionIndexMap.end()) return &sectionProperties[it->second];
        return nullptr;
    }

    const AreaInfo* getAreaInfo(int secId) const {
        auto it = areaIndexMap.find(secId);
        if (it != areaIndexMap.end()) return &areaInfos[it->second];
        return nullptr;
    }

    const SupportInfo* getSupportInfo(int nodeId) const {
        auto it = supportInfoMap.find(nodeId);
        if (it != supportInfoMap.end()) return &supportInfos[it->second];
        return nullptr;
    }

    const SpringSupport* getSpringSupport(int nodeId) const {
        for (const auto& sp : springSupports) {
            if (sp.nodeId == nodeId) return &sp;
        }
        return nullptr;
    }

    bool isSupportNode(int nodeId) const {
        if (!supportInfos.empty()) {
            return supportInfoMap.find(nodeId) != supportInfoMap.end();
        }
        for (const auto& sp : springSupports) {
            if (sp.nodeId == nodeId) return true;
        }
        return false;
    }

    double getMemberLength(const Element& elem) const {
        const Node* n1 = getNode(elem.n1);
        const Node* n2 = getNode(elem.n2);
        if (!n1 || !n2) return 0.0;
        double dx = n1->x - n2->x;
        double dy = n1->y - n2->y;
        double dz = n1->z - n2->z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    MemberType getMemberType(const Element& elem) const {
        const Node* n1 = getNode(elem.n1);
        const Node* n2 = getNode(elem.n2);
        if (!n1 || !n2) return MemberType::WebMember;
        if (n1->layer == 1 && n2->layer == 1) return MemberType::TopChord;
        if (n1->layer == 2 && n2->layer == 2) return MemberType::BottomChord;
        return MemberType::WebMember;
    }

    double getMemberWeight(const Element& elem) const {
        double L = getMemberLength(elem); // m
        const AreaInfo* area = getAreaInfo(elem.secId);
        if (!area) return 0.0;
        // 截面积 mm², 密度 7.85e-6 kg/mm³
        // 单位长度重量 (kg/m) = A (mm²) * 1e-6 (m²/mm²) * 7850 (kg/m³) = A * 7.85e-3 kg/m
        double linearMass = area->A * 7.85e-3;
        return L * linearMass; // kg
    }
};

} // namespace MST
