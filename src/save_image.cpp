#include <GLFW/glfw3.h>
#include <vector>
#include <algorithm>
#include <lodepng.h>
//
void save_image( GLFWwindow* window, std::string path ) {
	//
	int width, height;
	glfwGetFramebufferSize(window,&width,&height);
	std::vector<unsigned char> buffer(4*width*height);
	//
	glFlush();
	glPixelStorei(GL_UNPACK_ALIGNMENT,1);
	glReadPixels(0,0,width,height,GL_RGBA,GL_UNSIGNED_BYTE,buffer.data());
	//
	for(int j=0; j!=height/2; ++j) {
		std::swap_ranges(
			buffer.begin()+4*width*j,
			buffer.begin()+4*width*(j+1),
			buffer.begin()+4*width*(height-j-1));
	}
	//
	lodepng_encode32_file(path.c_str(),buffer.data(),width,height);
}