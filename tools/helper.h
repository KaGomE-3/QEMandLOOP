#pragma once
#include "tgaimage.h"
#include <iostream>
#include "model.h"
#include <vector>
#include <Eigen/Dense>
#include<Eigen/Core>

using namespace Eigen;
using namespace std;

Eigen::Matrix4f camera_matrix(Vec3f camera_point, Vec3f eye_point, Vec3f up) {
    // 
    Vec3f g = (camera_point-eye_point).normalize();
    Vec3f x = (up^g).normalize();
    Vec3f y = (g^x).normalize();
    Eigen::Matrix<float,4,4> r_view;
    // 世界坐标系为正坐标系，此旋旋转矩阵可以将（1，0，0）（世界）变到x
    r_view <<
        x.x, y.x, g.x, 0,
        x.y, y.y, g.y, 0,
        x.z, y.z, g.z, 0,
        0, 0, 0, 1;
    Eigen::Matrix4f r_view_inv = r_view.inverse();
    Eigen::Matrix<float,4,4> t_view;
    t_view <<
        1, 0, 0, -camera_point.x,
        0, 1, 0, -camera_point.y,
        0, 0, 1, -camera_point.z,
        0, 0, 0, 1;
    return r_view_inv * t_view;
}


Eigen::Matrix4f perspective_matrix() {
    float l = -0.5, r = 0.5, b = -0.5 ,t = 0.5, n = -0.5, f = -5.5;
    Eigen::Matrix<float, 4, 4> p1;
    p1 <<
        n, 0, 0, 0,
        0, n, 0, 0,
        0, 0, n + f, -n * f,
        0, 0, 1, 0;
    Eigen::Matrix<float, 4, 4> o1;
    o1 <<
        2/(r-l), 0, 0, 0,
        0, 2/(t-b), 0, 0,
        0, 0, 2 / (n - f), 0,
        0, 0,0, 1;
    Eigen::Matrix<float, 4, 4> o2;
    o2 <<
        1, 0, 0,-(r+l)/2,
        0, 1, 0, -(t+b)/2,
        0, 0, 1, -(n + f)/2,
        0, 0, 0, 1;
    Eigen::Matrix<float, 4, 4> tem = o1 * o2;
    Eigen::Matrix<float, 4, 4> res = tem*p1;

    return res;
}


Vec3f point2vp(Eigen::Matrix<float,4,4> v, Eigen::Matrix<float, 4, 4> p, Vec3f point) {
    Vector4f point_hom(point.x, point.y, point.z, 1.);
   
    point_hom = v * point_hom;
    /*cout << p << endl << endl << point_hom << endl << endl;*/
    point_hom = p * point_hom;
    //cout << point_hom << endl << endl;
    point_hom =  point_hom/point_hom[3];
    //cout << point_hom << endl << endl;
    return Vec3f(point_hom[0], point_hom[1], point_hom[2]);
}

float compute_w(Eigen::Matrix<float, 4, 4> v, Eigen::Matrix<float, 4, 4> p, Vec3f point) {
    Vector4f point_hom(point.x, point.y, point.z, 1.);

    point_hom = v * point_hom;
    point_hom = p * point_hom;
    
    return point_hom[3];
}

void line(Vec2f p1, Vec2f p2, TGAImage& image, TGAColor color) {
    // Bresenham算法，直线在x方向上每次增量为\Delta x=1，在y方向上每次的增量为\Delta y=k。
    //通过一个变量d将y方向上的累计增量记录下来，当d大于1时，标记点m进1，并对变量d进行-1操
    //作使得d的范围永远保持在[0,1]之间。并根据d的范围确定最终的y值，当d\leq 0.5时，则y保持不变；当d>0.5时，y加1
    bool exchange = 0;
    int x0 = p1.x, y0 = p1.y, x1 = p2.x, y1 = p2.y;
    // 坐标系xy互换，使k为-1 - 1，需记录，填充时需要坐标系互换
    if (abs(y1 - y0) > abs(x1 - x0)) {
        swap(x1, y1);
        swap(x0, y0);
        exchange = 1;
    }
    
    // 填充顺序不重要，无需记录
    if (x1 < x0) {
        swap(x1, x0);
        swap(y0, y1);
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
        d += abs(2*dy);
        // 因为像素点在中心而不是左下角顶点，所以超过0.5就上移，且包括中心点下部的y，因此要超过0.5之后+-1
        if (d > dx) {
            y1 > y0 ? y += 1 : y -= 1;
            d -= 2*dx;
        }
    }
}


Vec3f barycentric(Vec2f v1, Vec2f v2, Vec2f v3, Vec2f p)
{
    //别丢了分母等于0的情况
    if ((-(v1.x - v2.x) * (v3.y - v2.y) + (v1.y - v2.y) * (v3.x - v2.x)) == 0)
        return Vec3f(1, 0, 0);
    if (-(v2.x - v3.x) * (v1.y - v3.y) + (v2.y - v3.y) * (v1.x - v3.x) == 0)
        return Vec3f(1, 0, 0);
    float alpha = (-(p.x - v2.x) * (v3.y - v2.y) + (p.y - v2.y) * (v3.x - v2.x)) / (-(v1.x - v2.x) * (v3.y - v2.y) + (v1.y - v2.y) * (v3.x - v2.x));
    float beta = (-(p.x - v3.x) * (v1.y - v3.y) + (p.y - v3.y) * (v1.x - v3.x)) / (-(v2.x - v3.x) * (v1.y - v3.y) + (v2.y - v3.y) * (v1.x - v3.x));
    float gamma = 1 - alpha - beta;

    return Vec3f(alpha, beta, gamma);
}

bool IsInsideTriangle(Vec2f A, Vec2f B, Vec2f C, Vec2f p)
{
    double AB = (B.x - A.x) * (p.y - A.y) - (B.y - A.y) * (p.x - A.x);
    double BC = (C.x - B.x) * (p.y - B.y) - (C.y - B.y) * (p.x - B.x);
    double CA = (A.x - C.x) * (p.y - C.y) - (A.y - C.y) * (p.x - C.x);

    if (((AB > 0 && BC > 0 && CA > 0) || (AB < 0 && BC < 0 && CA < 0)))
        return true;
    else
        return false;
}

vector<vector<float>> msaa_check(Vec2f A, Vec2f B, Vec2f C, Vec2f p)
{
    vector<vector<float>>msaa_buffer;
    for (float i = -0.25; i <= 0.25; i += 0.5) {
        for (float j = -0.25; j <= 0.25; j += 0.5) {
            double AB = (B.x - A.x) * (p.y+i - A.y) - (B.y - A.y) * (p.x+j - A.x);
            double BC = (C.x - B.x) * (p.y+i - B.y) - (C.y - B.y) * (p.x+j - B.x);
            double CA = (A.x - C.x) * (p.y+i - C.y) - (A.y - C.y) * (p.x+j - C.x);

            if (((AB > 0 && BC > 0 && CA > 0) || (AB < 0 && BC < 0 && CA < 0))) {
                msaa_buffer.push_back({ 1, i, j });
            }
            else {
                msaa_buffer.push_back({ 0, i, j });
            }
        }
    }
    return msaa_buffer;
}



// 光栅化、着色函数
void triangle_render(Model* model,vector<Vec3f>ps, vector<Vec2i>uvs, TGAImage& image, vector<float>&buffer,
    int weith,int height, vector<Vec3f>ns, Vec3f camera_point,Vec3f light_dir,
    vector<Vec3f>ps_prior, vector<float>ws) {

    // 环境光
    float ambient = 0.1;
    
    // 到图像坐标
    float x1 = (ps[0].x + 1.) * weith / 2., y1 = (ps[0].y + 1.) * height / 2.;
    float x2 = (ps[1].x + 1.) * weith / 2. , y2 = (ps[1].y + 1.) * height / 2.;
    float x3 = (ps[2].x + 1.) * weith / 2. , y3 = (ps[2].y + 1.) * height / 2.;

    // 图像坐标系包围框
    float x_min = min(x1, x2);
    x_min = min(x_min, x3);

    float x_max = max(x1, x2);
    x_max = max(x_max, x3);

    float y_min = min(y1, y2);
    y_min = min(y_min, y3);

    float y_max = max(y1, y2);
    y_max = max(y_max, y3);

    for (int x = x_min-1; x <= x_max+1; x++) {
        for (int y = y_min-1; y <= y_max+1; y++) {
            if (x < 0 || x > weith || y < 0 || y > height) {
                continue;
            }
            if (IsInsideTriangle(Vec2f(x1, y1), Vec2f(x2, y2), Vec2f(x3, y3), Vec2f(x + 0.5, y + 0.5))) {
                Vec3f bary_centric = barycentric(Vec2f(x1, y1), Vec2f(x2, y2), Vec2f(x3, y3), Vec2f(x + 0.5, y + 0.5));
                // 透视矫正
                float k = 1.0 / (bary_centric.x * ws[0] + bary_centric.y * ws[1] + bary_centric.z * ws[2]);
                float alpha = bary_centric.x / (ws[0] * k);
                float beta = bary_centric.y / (ws[1] * k);
                float gama = bary_centric.z / (ws[2] * k);

                float p_depth = (alpha * ps[0].z + beta * ps[1].z + gama * ps[2].z) / 2.0 + 0.5;
                int idx = y * height + x;
                if (buffer[idx] < p_depth) {
                    buffer[idx] = p_depth;
                    Vec3f interpertion_point = ps_prior[0] * alpha + ps_prior[1] * beta + ps_prior[2] * gama;

                    Vec3f eye_light = camera_point - interpertion_point;
                    eye_light.normalize();
                    Vec3f normal = ns[0] * alpha + ns[1] * beta + ns[2] * gama;
                    normal.normalize();
                    light_dir.normalize();
                    Vec3f half_v = (eye_light + light_dir);
                    half_v.normalize();

                    float blin_phong = pow(max(0.f, half_v * normal), 16);
                    float indensity = max(0.f, normal * light_dir);
                    Vec2i color_index = Vec2i(alpha * uvs[0].x + beta * uvs[1].x + gama * uvs[2].x, alpha * uvs[0].y + beta * uvs[1].y + gama * uvs[2].y);
                    TGAColor color = model->diffuse(color_index);

                    unsigned char r = min(int(color.r * (indensity + ambient + blin_phong)), 255);
                    unsigned char g = min(int(color.g * (indensity + ambient + blin_phong)), 255);
                    unsigned char b = min(int(color.b * (indensity + ambient + blin_phong)), 255);

                    TGAColor res_color = TGAColor(r, g, b, 255);
                    image.set(x, y, res_color);
                }
            }
        }
    }
}

//画非洲老哥的头
void drawAfrican(Model* model, TGAImage& image,Vec3f light, int weith, int height,
    Eigen::Matrix4f camera_matrix, Eigen::Matrix4f p_matrix, Vec3f camera_point) {
    // 深度缓存
    vector<float>buffer(height* weith, -99999999999.);
    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i); 
        Vec3f p1_prior = model->vert(face[0]);
        Vec3f p2_prior = model->vert(face[1]);
        Vec3f p3_prior = model->vert(face[2]);
        Vec2i uv1(0, 0);
        Vec2i uv2(0, 0);
        Vec2i uv3(0, 0);
        uv1 = model->uv(i, 0);
        uv2 = model->uv(i, 1);
        uv3 = model->uv(i, 2);
        Vec3f n1 = model->normal(i, 0);
        Vec3f n2 = model->normal(i, 1);
        Vec3f n3 = model->normal(i, 2);
        light.normalize();
        Vec3f p1 = point2vp(camera_matrix, p_matrix,  p1_prior);
        Vec3f p2 = point2vp(camera_matrix, p_matrix,p2_prior);
        Vec3f p3 = point2vp(camera_matrix, p_matrix, p3_prior);

        float w1 = compute_w(camera_matrix, p_matrix, p1_prior);
        float w2 = compute_w(camera_matrix, p_matrix, p2_prior);
        float w3 = compute_w(camera_matrix, p_matrix, p3_prior);

        vector<Vec3f>ps{ p1,p2,p3 };
        vector<Vec2i>uvs{ uv1,uv2,uv3 };
        vector<float>ws{w1,w2,w3 };
        vector<Vec3f>ps_prior{p1_prior,p2_prior,p3_prior};
        vector<Vec3f>ns{n1, n2, n3};

        triangle_render(model, ps, uvs, image,  buffer,weith,height,ns, camera_point,light,ps_prior,ws);

    }
}


