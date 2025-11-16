#include "tools/tgaimage.h"
#include "tools/model.h"

#include <vector>
#include <map>
#include <utility>
#include <algorithm>  
#include <cmath>


const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);

// Bresenham 画线函数（保持你原来的）
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
        d += std::abs(2 * dy);
        if (d > dx) {
            y1 > y0 ? y += 1 : y -= 1;
            d -= 2 * dx;
        }
    }
}

// 简单三角面结构
struct Face {
    int v[3];
};

// 边信息：记录对边顶点和新生成的边点索引
struct EdgeInfo {
    int vOpp[2];   // 两个对顶点（内部边）
    int newIndex;  // 新生成的边点顶点索引
    EdgeInfo() {
        vOpp[0] = vOpp[1] = -1;
        newIndex = -1;
    }
};

static std::pair<int, int> make_edge(int a, int b) {
    if (a > b) std::swap(a, b);
    return std::make_pair(a, b);
}

// 单次 Loop 细分
void loopSubdivisionOnce(std::vector<Vec3f>& verts, std::vector<Face>& faces) {
    int nVerts = (int)verts.size();
    int nFaces = (int)faces.size();

    // 1. 顶点邻接（用于更新老顶点）
    std::vector<std::vector<int>> vertexNeighbors(nVerts);
    // 2. 边信息（用于生成边点）
    std::map<std::pair<int, int>, EdgeInfo> edges;

    // 构建邻接表和边的对顶点信息
    for (int i = 0; i < nFaces; ++i) {
        int i0 = faces[i].v[0];
        int i1 = faces[i].v[1];
        int i2 = faces[i].v[2];

        // 邻接
        auto addNeighbor = [&](int a, int b) {
            vertexNeighbors[a].push_back(b);
            };
        addNeighbor(i0, i1); addNeighbor(i1, i0);
        addNeighbor(i1, i2); addNeighbor(i2, i1);
        addNeighbor(i2, i0); addNeighbor(i0, i2);

        // 每条边记录对顶点
        auto processEdge = [&](int a, int b, int opp) {
            auto key = make_edge(a, b);
            EdgeInfo& info = edges[key];
            if (info.vOpp[0] == -1) info.vOpp[0] = opp;
            else if (info.vOpp[1] == -1) info.vOpp[1] = opp;
            };

        processEdge(i0, i1, i2);
        processEdge(i1, i2, i0);
        processEdge(i2, i0, i1);
    }

    // 去重邻接（防止重复）
    for (int i = 0; i < nVerts; ++i) {
        auto& nbrs = vertexNeighbors[i];
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    // 3. 计算旧顶点的新位置
    std::vector<Vec3f> newVerts;
    newVerts.resize(nVerts); // 先只放老顶点的新位置

    for (int i = 0; i < nVerts; ++i) {
        std::vector<int>& nbrs = vertexNeighbors[i];
        int n = (int)nbrs.size();
        if (n == 0) { // 孤立点（理论上不会出现）
            newVerts[i] = verts[i];
            continue;
        }

        float beta;
        if (n == 3) beta = 3.0f / 16.0f;
        else        beta = 3.0f / (8.0f * n);

        Vec3f sumNeighbors(0, 0, 0);
        for (int j = 0; j < n; ++j) {
            sumNeighbors = sumNeighbors + verts[nbrs[j]];
        }
        float scaleSelf = 1.0f - n * beta;
        newVerts[i] = verts[i] * scaleSelf + sumNeighbors * beta;
    }

    // 4. 为每条边生成新的顶点（边点）
    for (auto& it : edges) {
        const std::pair<int, int>& e = it.first;
        EdgeInfo& info = it.second;
        int a = e.first;
        int b = e.second;
        Vec3f v0 = verts[a];
        Vec3f v1 = verts[b];

        Vec3f newPos;
        if (info.vOpp[0] != -1 && info.vOpp[1] != -1) {
            // 内部边：3/8 * (v0 + v1) + 1/8 * (v2 + v3)
            Vec3f v2 = verts[info.vOpp[0]];
            Vec3f v3 = verts[info.vOpp[1]];
            newPos = (v0 + v1) * (3.0f / 8.0f) + (v2 + v3) * (1.0f / 8.0f);
        }
        else {
            // 边界边（如果有）：简单取中点
            newPos = (v0 + v1) * 0.5f;
        }

        info.newIndex = (int)newVerts.size();
        newVerts.push_back(newPos);
    }

    // 5. 生成新的三角面（每个旧三角拆成 4 个）
    std::vector<Face> newFaces;
    newFaces.reserve(nFaces * 4);

    for (int i = 0; i < nFaces; ++i) {
        int i0 = faces[i].v[0];
        int i1 = faces[i].v[1];
        int i2 = faces[i].v[2];

        int e01 = edges[make_edge(i0, i1)].newIndex;
        int e12 = edges[make_edge(i1, i2)].newIndex;
        int e20 = edges[make_edge(i2, i0)].newIndex;

        Face f1{ { i0, e01, e20 } };
        Face f2{ { i1, e12, e01 } };
        Face f3{ { i2, e20, e12 } };
        Face f4{ { e01, e12, e20 } };

        newFaces.push_back(f1);
        newFaces.push_back(f2);
        newFaces.push_back(f3);
        newFaces.push_back(f4);
    }

    // 替换原来的顶点和面
    verts.swap(newVerts);
    faces.swap(newFaces);
}

int main(int argc, char** argv) {
    Model model{ "D:/graphics/test/test/asset/african_head.obj" };

    // 取出模型原始网格
    std::vector<Vec3f> verts;
    verts.reserve(model.nverts());
    for (int i = 0; i < model.nverts(); ++i) {
        verts.push_back(model.vert(i));
    }

    std::vector<Face> faces;
    faces.reserve(model.nfaces());
    for (int i = 0; i < model.nfaces(); ++i) {
        std::vector<int> face = model.face(i); // 这里假设 face.size() == 3
        Face f{ { face[0], face[1], face[2] } };
        faces.push_back(f);
    }

    // 执行一次 Loop 细分（可以改成多次）
    int iterations = 1;
    for (int it = 0; it < iterations; ++it) {
        loopSubdivisionOnce(verts, faces);
    }

    // 渲染细分后的网格线框
    int width = 800;
    int height = 800;
    TGAImage image(width, height, TGAImage::RGB);

    for (int i = 0; i < (int)faces.size(); i++) {
        for (int j = 0; j < 3; j++) {
            int idx0 = faces[i].v[j];
            int idx1 = faces[i].v[(j + 1) % 3];

            Vec3f v0 = verts[idx0];
            Vec3f v1 = verts[idx1];

            int x0 = (int)((v0.x + 1.f) * width / 2.f);
            int y0 = (int)((v0.y + 1.f) * height / 2.f);
            int x1 = (int)((v1.x + 1.f) * width / 2.f);
            int y1 = (int)((v1.y + 1.f) * height / 2.f);

            line(x0, y0, x1, y1, image, white);
        }
    }

    image.flip_vertically();
    image.write_tga_file("complexoutput.tga");
    return 0;
}
