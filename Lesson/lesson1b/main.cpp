#include <cmath>
#include "tgaimage.h"
#include <vector>
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(x0-x1)<std::abs(y0-y1)){
        std::swap(x0,y0);
        std::swap(x1,y1);
        steep=true;
    }
    if(x0>x1){
        std::swap(x0,x1);
        std::swap(y0,y1);
        
    }
    int dx = x1-x0;
    int dy = y1-y0;
    int derror2 = std::abs(dy)*2;
    int error2 = 0;
//    float derror = std::abs(dy/float(dx)); //use float, optimal for graphics
//    float error = 0;
    int y = y0;
    for (int x = x0; x<= x1; x++){
        if (steep){
            image.set(y,x,color);
        } else {
            image.set(x,y,color);
        }
        error2 +=derror2;
        if (error2>dx){
            y+=(y1>y0?1:-1);
            error2 -= dx*2;
        }
    }

}


int main(int argc, char** argv) {
	TGAImage image(height*2, width*2, TGAImage::RGB);
	//image.set(52, 41, red);// adds a red dot
	//line(13,20,80,40,image,white);
//    for (int i = 0; i<10; i++){
//        line(13,20,80,40,image,white);
//        line(20,13,40,40,image, red);
//        line(80,40,13,20,image,red);
//    }
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }
//    TGAImage image(width, height, TGAImage::RGB);
    for (int i=0; i<model->nfaces();i++){
        std::vector<int> face = model->face(i);
        for (int j=0; j<3; j++){
            Vec3f v0 = model->vert(face[j]);
            Vec3f v1 = model->vert(face[(j+1)%3]);
            int x0 = (v0.x+1.)*width/2.;
            int y0 = (v0.y+1.)*height/2.;
            int x1 = (v1.x+1.)*width/2.;
            int y1 = (v1.y+1.)*height/2.;
            line(x0,y0,x1,y1,image,white);
        }
    
    }
    
	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
    delete model;
	return 0;
}

