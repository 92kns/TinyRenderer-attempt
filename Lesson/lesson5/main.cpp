#include <cmath>
#include "tgaimage.h"
#include <vector>
#include "model.h"
#include <cstdlib>
#include "geometry.h"
#include <limits>
#include <iostream>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const TGAColor blue  = TGAColor(0,   0,   255, 255);
Model *model = NULL;
const int width  = 800;
const int height = 800;
const int depth = 255;
int *zbuffer = NULL;
Vec3f light_dir = Vec3f(1,-1,1).normalize();
Vec3f eye(1,1,3); // test coordinates
Vec3f center(0,0,0);

Matrix viewport(int x, int y, int w, int h){
    Matrix m = Matrix::identity(4);
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;
    
    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}

Matrix lookat( Vec3f eye, Vec3f center, Vec3f up){
    Vec3f z = (eye-center).normalize();
    Vec3f x = (up^z).normalize();
    Vec3f y = (z^x).normalize();
    Matrix res = Matrix::identity(4);
    for (int i = 0; i < 3 ; i++){
        res[0][i] = x[i];
        res[1][i] = y[i];
        res[2][i] = z[i];
        res[i][3] = -center[i];
    }
    return res;
}

void line2(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(p0.x-p1.x)<std::abs(p0.y-p1.y)) {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }
    if (p0.x>p1.x) {
        std::swap(p0, p1);
    }

    for (int x=p0.x; x<=p1.x; x++) {
        float t = (x-p0.x)/(float)(p1.x-p0.x);
        int y = p0.y*(1.-t) + p1.y*t+.5;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color){
    bool steep = false;
    if (std::abs(x0-x1) < std::abs(y0-y1)){
        std::swap(x0,y0);
        std::swap(x1,y1);
        steep = true;
    }
    if(x0>x1) {
        std::swap(x0,x1);
        std::swap(y0,y1);
        
    }
    
    for (int x = x0; x<= x1; x++){
        float t= (x-x0)/(float)(x1-x0);
        int y = y0*(1.-t) + y1*t;
        if (steep) {
            image.set(y,x,color);
        } else {
            image.set(x,y,color);
        }
    }
}

//Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P){
//    Vec3f s[2];
//    for ( int i = 2 ; i --;){
//        s[i][0] = C[i]-A[i];
//        s[i][1] = B[i]-A[i];
//        s[i][2] = A[i]-P[i];
//    }
//    Vec3f u = cross(s[0], s[1]);
//    if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
//        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
//    return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
//    
//}

//void triangleold(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color){
//    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
//    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
//    Vec2f clamp(image.get_width()-1, image.get_height()-1);
//    
//    for (int i =0; i <3 ; i++){
//        for(int j = 0; j<2; j++){
//            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
//            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
//        }
//    }
//    
//    Vec3f P;
//    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++){
//        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++){
//            Vec3f bc_screen = barycentric(pts[0], pts[1], pts[2], P);
//            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
//            P.z = 0;
//            for ( int i =0; i <3; i++) P.z += pts[i][2]*bc_screen[i];
//            if(zbuffer[int(P.x+P.y*width)]<P.z){
//                zbuffer[int(P.x+P.y*width)] = P.z;
//                image.set(P.x, P.y, color);
//        
//        }
//        }
//    }
//}


void triangle(Vec3i t0, Vec3i t1, Vec3i t2, float ity0, float ity1, float ity2, TGAImage &image, int *zbuffer) {
    if (t0.y==t1.y && t0.y==t2.y) return; // i dont care about degenerate triangles
    if (t0.y>t1.y) { std::swap(t0, t1); std::swap(ity0, ity1); }
    if (t0.y>t2.y) { std::swap(t0, t2); std::swap(ity0, ity2); }
    if (t1.y>t2.y) { std::swap(t1, t2); std::swap(ity1, ity2); }
    
    int total_height = t2.y-t0.y;
    for (int i=0; i<total_height; i++) {
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here
        Vec3i A    =               t0  + Vec3f(t2-t0  )*alpha;
        Vec3i B    = second_half ? t1  + Vec3f(t2-t1  )*beta : t0  + Vec3f(t1-t0  )*beta;
        float ityA =               ity0 +   (ity2-ity0)*alpha;
        float ityB = second_half ? ity1 +   (ity2-ity1)*beta : ity0 +   (ity1-ity0)*beta;
        if (A.x>B.x) { std::swap(A, B); std::swap(ityA, ityB); }
        for (int j=A.x; j<=B.x; j++) {
            float phi = B.x==A.x ? 1. : (float)(j-A.x)/(B.x-A.x);
            Vec3i    P = Vec3f(A) +  Vec3f(B-A)*phi;
            float ityP =    ityA  + (ityB-ityA)*phi;
            int idx = P.x+P.y*width;
            if (P.x>=width||P.y>=height||P.x<0||P.y<0) continue;
            if (zbuffer[idx]<P.z) {
                zbuffer[idx] = P.z;
                image.set(P.x, P.y, TGAColor(255, 255, 255)*ityP);
            }
        }
    }
}

void triangle_mine(Vec3i t0, Vec3i t1, Vec3i t2, float ity0, float ity1, float ity2, TGAImage &image, int *zbuffer) {
    // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!)
    if (t0.y>t1.y) { std::swap(t0, t1); std::swap(ity0, ity1); }
    if (t0.y>t2.y) { std::swap(t0, t2); std::swap(ity0, ity2); }
    if (t1.y>t2.y) { std::swap(t1, t2); std::swap(ity1, ity2); }
    
    

    int total_height = t2.y-t0.y;
    for (int y=t0.y; y<=t2.y; y++) {

        bool firsthalf = y<t1.y; // || t0.y == t1.y;
        int segment_height = firsthalf? t1.y-t0.y: t2.y - t1.y;
        float alpha = (float)(y-t0.y)/total_height;
        float beta  = (float)(y- (firsthalf? t0.y:t1.y ))/segment_height; // be careful with divisions by zero
        Vec3i A = t0 + Vec3f(t2-t0)*alpha;
        Vec3i B = firsthalf? t0 + Vec3f(t1-t0)*beta : t1+ Vec3f(t2-t1)*beta;
        
        float ityA = ity0 + (ity2-ity0)*alpha;
        float ityB = firsthalf?ity0 + (ity1-ity0)*beta :ity1 + (ity2-ity1)*beta;
        
        if (A.x>B.x) {std::swap(A, B); std::swap(ityA,ityB);}
//        std::cout<<"made it here"<<std::endl;
        for (int j=A.x; j<=B.x; j++) {
//            image.set(j, y, color); // attention, due to int casts t0.y+i != A.y
            float phi = B.x == A.x ? 1. : (float)(j-A.x)/(B.x - A.x);
            Vec3i P = Vec3f(A) + Vec3f(B-A)*phi;
            float ityP = ityA + (ityB - ityA)*phi;
            int idx = P.x + P.y*width;
//            std::cout<<"made it here"<<std::endl;
            if(P.x>=width||P.y>=height||P.x<0||P.y<0) continue;
            if(zbuffer[idx]<P.z){
                zbuffer[idx] = P.z;
                image.set(P.x, P.y, TGAColor(255,255,255)*ityP);
//                std::cout<<"made it here"<<std::endl;
            }
        }
    }
    
}

Vec3f world2screen(Vec3f v){
    return Vec3f(int((v.x+1.)*width/2.+.5), int((v.y+1.)*height/2.+.5),v.z);
}

//void triangle2(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
//    line(t0, t1, image, color);
//    line(t1, t2, image, color);
//    line(t2, t0, image, color);
//}

void triangle3(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
    // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!)
    if (t0.y>t1.y) std::swap(t0, t1);
    if (t0.y>t2.y) std::swap(t0, t2);
    if (t1.y>t2.y) std::swap(t1, t2);
    int total_height = t2.y-t0.y;
    for (int y=t0.y; y<=t1.y; y++) {
        int segment_height = t1.y-t0.y+1;
        float alpha = (float)(y-t0.y)/total_height;
        float beta  = (float)(y-t0.y)/segment_height; // be careful with divisions by zero
        Vec2i A = t0 + (t2-t0)*alpha;
        Vec2i B = t0 + (t1-t0)*beta;
        if (A.x>B.x) std::swap(A, B);
        for (int j=A.x; j<=B.x; j++) {
            image.set(j, y, color); // attention, due to int casts t0.y+i != A.y
        }
    }
    for (int y=t1.y; y<=t2.y; y++) {
        int segment_height =  t2.y-t1.y+1;
        float alpha = (float)(y-t0.y)/total_height;
        float beta  = (float)(y-t1.y)/segment_height; // be careful with divisions by zero
        Vec2i A = t0 + (t2-t0)*alpha;
        Vec2i B = t1 + (t2-t1)*beta;
        if (A.x>B.x) std::swap(A, B);
        for (int j=A.x; j<=B.x; j++) {
            image.set(j, y, color); // attention, due to int casts t0.y+i != A.y
        }
    }
}

void triangle4(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
    // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!)
    if (t0.y>t1.y) std::swap(t0, t1);
    if (t0.y>t2.y) std::swap(t0, t2);
    if (t1.y>t2.y) std::swap(t1, t2);
    int total_height = t2.y-t0.y;
    for (int y=t0.y; y<=t2.y; y++) {
        bool firsthalf = y<t1.y; // || t0.y == t1.y;
        int segment_height = firsthalf? t1.y-t0.y: t2.y - t1.y;
        float alpha = (float)(y-t0.y)/total_height;
        float beta  = (float)(y- (firsthalf? t0.y:t1.y ))/segment_height; // be careful with divisions by zero
        Vec2i A = t0 + (t2-t0)*alpha;
        Vec2i B = firsthalf? t0 + (t1-t0)*beta: t1+(t2-t1)*beta;
        if (A.x>B.x) std::swap(A, B);
        for (int j=A.x; j<=B.x; j++) {
            image.set(j, y, color); // attention, due to int casts t0.y+i != A.y
        }
    }
//    for (int y=t1.y; y<=t2.y; y++) {
//        int segment_height =  t2.y-t1.y+1;
//        float alpha = (float)(y-t0.y)/total_height;
//        float beta  = (float)(y-t1.y)/segment_height; // be careful with divisions by zero
//        Vec2i A = t0 + (t2-t0)*alpha;
//        Vec2i B = t1 + (t2-t1)*beta;
//        if (A.x>B.x) std::swap(A, B);
//        for (int j=A.x; j<=B.x; j++) {
//            image.set(j, y, color); // attention, due to int casts t0.y+i != A.y
//        }
//    }
}


void rasterize(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color, int ybuffer[])
{
    if (p0.x>p1.x){
        std::swap(p0,p1);
    }
    for (int x=p0.x;x<= p1.x; x++){
        float t = (x-p0.x)/(float)(p1.x-p0.x);
        int y = p0.y*(1.-t) +p1.y*t + .5;
        if (ybuffer[x]<y) {
            ybuffer[x] = y;
            image.set(x,0,color);
        }
    }
}
int main(int argc, char** argv) {
//    TGAImage image(width, height, TGAImage::RGB);

//    Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)};
//    Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)};
//    Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)};
//
//    triangle(t0[0], t0[1], t0[2], image, red);
//    triangle(t1[0], t1[1], t1[2], image, white);
//    triangle(t2[0], t2[1], t2[2], image, green);
//    { // just dumping the 2d scene (yay we have enough dimensions!)
//        TGAImage scene(width, height, TGAImage::RGB);
//        
//        // scene "2d mesh"
//        line(Vec2i(20, 34),   Vec2i(744, 400), scene, red);
//        line(Vec2i(120, 434), Vec2i(444, 400), scene, green);
//        line(Vec2i(330, 463), Vec2i(594, 200), scene, blue);
//        
//        // screen line
//        line(Vec2i(10, 10), Vec2i(790, 10), scene, white);
//        
//        scene.flip_vertically(); // i want to have the origin at the left bottom corner of the image
//        scene.write_tga_file("scene.tga");
//    }
//    
//    {
//        TGAImage render(width, 16, TGAImage::RGB);
//        int ybuffer[width];
//        for (int i = 0; i <width; i++){
//            ybuffer[i] = std::numeric_limits<int>::min();
//        }
//        rasterize(Vec2i(20, 34),   Vec2i(744, 400), render, red,   ybuffer);
//        rasterize(Vec2i(120, 434), Vec2i(444, 400), render, green, ybuffer);
//        rasterize(Vec2i(330, 463), Vec2i(594, 200), render, blue,  ybuffer);
//        
//        for ( int i = 0; i <width; i++){
//            for(int j = 1; j < 16; j++){
//                render.set(i,j,render.get(i,0));
//                
//            }
//        }
//        
//        render.flip_vertically(); // i want to have the origin at the left bottom corner of the image
//        render.write_tga_file("render.tga");
//    }
    
    if(2 == argc){
        model = new Model(argv[1]);
        
    } else {
        model = new Model("obj/african_head.obj");
    }

    zbuffer = new int[width*height];
    for ( int i =0; i< width*height; i++){
        zbuffer[i] = std::numeric_limits<float>::min();
    }
    
  //draw
    {
        Matrix ModelView = lookat(eye, center, Vec3f(0,1,0));
        Matrix Projection = Matrix::identity(4);
        Matrix ViewPort = viewport(width/8, height/8, width*3/4, height*3/4);
        Projection[3][2] = -1.f/(eye-center).norm();
        
        std::cerr << ModelView << std::endl;
        std::cerr << Projection << std::endl;
        std::cerr << ViewPort << std::endl;
        Matrix z  = (ViewPort*Projection*ModelView);
        std::cerr << z << std::endl;
        
        TGAImage image(width, height, TGAImage::RGB);
        for ( int i=0; i<model->nfaces(); i++){
            std::vector<int> face = model->face(i);
            Vec3i screen_coords[3];
            Vec3f world_coords[3];
            float intensity[3];
            for (int j = 0; j <3; j++){
                Vec3f v = model->vert(face[j]);
                screen_coords[j] = Vec3f(ViewPort*Projection*ModelView*Matrix(v));
                world_coords[j] = v;
                intensity[j] = model->norm(i,j)*light_dir;
            }
//            std::cout<<"made it here"<<std::endl;
           triangle(screen_coords[0], screen_coords[1], screen_coords[2], intensity[0], intensity[1], intensity[2], image, zbuffer);
           // triangle(screen_coords[0], screen_coords[1], screen_coords[2], intensity[0])
        }
//        std::cout<<"made it here"<<std::endl;
        image.flip_vertically();
        image.write_tga_file("outputlesson5.tga");
        
    }
    { // dump z-buffer (debugging purposes only)
        TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                zbimage.set(i, j, TGAColor(zbuffer[i+j*width]));
            }
        }
        zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
        zbimage.write_tga_file("zbuffer.tga");
    }
//    std::cout<<"made it here"<<std::endl;
    delete model;
    delete [] zbuffer;
    return 0;
}
