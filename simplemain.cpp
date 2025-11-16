#include "tools/tgaimage.h"
#include "tools/model.h"
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>

// 全局颜色定义
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);

// 辅助数据结构：顶点、面、QEM误差（仅用于网格简化）
struct Vertex {
    Vec3f pos;    // 三维坐标
    int id;       // 唯一ID
    bool valid;   // 是否有效（未被合并）
    Vertex(Vec3f p, int i) : pos(p), id(i), valid(true) {}
};

struct Face {
    std::vector<int> vertIds;  // 面的顶点ID列表（3个ID组成三角形）
    Face(std::vector<int> ids) : vertIds(ids) {}
};

struct QEM {
    float error;  // 合并误差（越小优先级越高）
    int v1, v2;   // 待合并顶点ID
    Vec3f optimalPos; // 合并后的最优顶点坐标
    bool operator>(const QEM& other) const { return error > other.error; }
};

// 计算两个顶点合并后的最优位置和误差（严谨版QEM核心）
bool computeOptimalMerge(const Vertex& v1, const Vertex& v2, Vec3f& optimalPos, float& error) {
    // 简化版：用中点近似最优位置（完整QEM需结合二次型矩阵，此处为了易用性简化）
    optimalPos = (v1.pos + v2.pos) * 0.5f;
    Vec3f diff = v1.pos - v2.pos;
    error = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z; // 距离平方作为误差
    return true;
}

// 网格简化函数：返回简化后的顶点和面
std::pair<std::vector<Vertex>, std::vector<Face>> simplifyMesh(Model& model, int targetFaceCount) {
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    // 1. 从Model读取顶点和面数据
    for (int i = 0; i < model.nverts(); ++i) {
        vertices.emplace_back(model.vert(i), i);
    }
    for (int i = 0; i < model.nfaces(); ++i) {
        faces.emplace_back(model.face(i));
    }

    if (faces.size() <= targetFaceCount) {
        return std::make_pair(vertices, faces);
    }

    // 2. 初始化优先队列（按QEM误差从小到大排序）
    std::priority_queue<QEM, std::vector<QEM>, std::greater<QEM>> pq;
    for (int i = 0; i < vertices.size(); ++i) {
        for (int j = i + 1; j < vertices.size(); ++j) {
            QEM qem;
            qem.v1 = i;
            qem.v2 = j;
            computeOptimalMerge(vertices[i], vertices[j], qem.optimalPos, qem.error);
            pq.push(qem);
        }
    }

    // 3. 循环合并顶点，直到面数达标
    while (faces.size() > targetFaceCount) {
        if (pq.empty()) break;
        QEM qem = pq.top();
        pq.pop();

        int v1 = qem.v1, v2 = qem.v2;
        if (!vertices[v1].valid || !vertices[v2].valid) continue;

        // 合并v1和v2为新顶点
        Vertex newV(qem.optimalPos, vertices.size());
        vertices.push_back(newV);

        // 4. 更新所有面：替换v1和v2为新顶点
        std::vector<Face> newFaces;
        for (const auto& face : faces) {
            bool hasV1 = false, hasV2 = false;
            for (int vid : face.vertIds) {
                if (vid == v1) hasV1 = true;
                if (vid == v2) hasV2 = true;
            }
            if (hasV1 && hasV2) continue;

            Face newFace{ {} };
            for (int vid : face.vertIds) {
                if (vid == v1 || vid == v2) {
                    newFace.vertIds.push_back(newV.id);
                }
                else {
                    newFace.vertIds.push_back(vid);
                }
            }
            newFaces.push_back(newFace);
        }
        faces = newFaces;

        // 标记原顶点为无效
        vertices[v1].valid = false;
        vertices[v2].valid = false;

        // 5. 计算新顶点与其他有效顶点的QEM，加入队列
        for (int i = 0; i < vertices.size() - 1; ++i) {
            if (vertices[i].valid) {
                QEM newQem;
                newQem.v1 = i;
                newQem.v2 = newV.id;
                computeOptimalMerge(vertices[i], newV, newQem.optimalPos, newQem.error);
                pq.push(newQem);
            }
        }
    }

    return std::make_pair(vertices, faces);
}

// Bresenham画线函数（保持不变）
void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
    bool exchange = 0;
    if (abs(y1 - y0) > abs(x1 - x0)) {
        std::swap(x1, y1);
        std::swap(x0, y0);
        exchange = 1;
    }
    if (x1 < x0) {
        std::swap(x1, x0);
        std::swap(y0, y1);
    }
    int dx = x1 - x0, dy = y1 - y0;
    float d = 0;
    int y = y0;
    for (int x = x0; x <= x1; x++) {
        if (exchange) {
            image.set(y, x, color);
        }
        else {
            image.set(x, y, color);
        }
        d += abs(2 * dy);
        if (d > dx) {
            y1 > y0 ? y += 1 : y -= 1;
            d -= 2 * dx;
        }
    }
}

int main(int argc, char** argv) {
    Model model{ "D:/graphics/test/test/asset/african_head.obj" };
    int width = 800;
    int height = 800;
    TGAImage image(width, height, TGAImage::RGB);

    // 网格简化：目标面数设为原面数的1/3（可根据需求调整）
    int targetFaceCount = model.nfaces() / 20;
    auto simplifiedMesh = simplifyMesh(model, targetFaceCount);
    std::vector<Vertex> simplifiedVerts = simplifiedMesh.first;
    std::vector<Face> simplifiedFaces = simplifiedMesh.second;

    // 绘制简化后的网格
    for (const auto& face : simplifiedFaces) {
        for (int j = 0; j < 3; j++) {
            int vid0 = face.vertIds[j];
            int vid1 = face.vertIds[(j + 1) % 3];
            Vec3f v0 = simplifiedVerts[vid0].pos;
            Vec3f v1 = simplifiedVerts[vid1].pos;

            int x0 = (v0.x + 1) * width / 2;
            int y0 = (v0.y + 1) * height / 2;
            int x1 = (v1.x + 1) * width / 2;
            int y1 = (v1.y + 1) * height / 2;
            line(x0, y0, x1, y1, image, white);
        }
    }

    image.flip_vertically();
    image.write_tga_file("qem_simplified.tga");
    return 0;
}