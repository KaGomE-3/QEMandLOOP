#include "tools/tgaimage.h"
#include "tools/model.h"
#include"vector"


const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
    bool exchange = 0;
    if (abs(y1 - y0) > abs(x1 - x0)) {
        std::swap(x1, y1);
        std::swap(x0, y0);
        exchange = 1;
    }
    // 填充顺序不重要，无需记录
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
        // 相关计算变量都*2dx，整数运算提升效率
        d += abs(2 * dy);
        // 因为像素点在中心而不是左下角顶点，所以超过0.5就上移，且包括中心点下部的y，因此要超过0.5之后+-1
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

    for (int i = 0; i < model.nfaces(); i++) {
        std::vector<int> face = model.face(i);
        for (int j = 0; j < 3; j++) {
            Vec3f v0 = model.vert(face[j]);
            Vec3f v1 = model.vert(face[(j + 1) % 3]);
            int x0 = (v0.x + 1) * width / 2;
            int y0 = (v0.y + 1) * height / 2;
            int x1 = (v1.x + 1) * width / 2;
            int y1 = (v1.y + 1) * height / 2;
            line(x0, y0, x1, y1, image, white);
        }
    }
    //image.set(52, 41, white);
    //line(2, 3, 52, 53, image, white);
    image.flip_vertically();
    image.write_tga_file("output.tga");

    return 0;
}