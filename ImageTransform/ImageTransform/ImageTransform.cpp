// ImageTransform.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "lodepng.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define THETA_MAX 2.0*M_PI
#define THETA_MIN 0.0
#define PHI_MAX 1.57
#define PHI_MIN -1.0487

int writeImage(int width, int height, vector<unsigned char>& buf, vector<double>& xyz = vector<double>(), vector<double>& depth = vector<double>(), double minDist = 0.5, double maxDist = 5);
int getIndex(Vector3d& point, int width, int height, vector<double>& blend, bool debug = false);
int getIndex(int width, int height, double x, double y);

int main(int argc, char* argv[])
{
	if (argc < 1) {
		cout << "Usage: ImageTransform infile\n  infile: the input file (*.ptx or *.png)" << endl;
		exit(1);
	}

	string infile = argv[1];

	if (infile.substr(infile.find_last_of('.')+1) == "ptx") {
		ifstream ifs;
		ifs.open(infile);

		if (ifs.is_open()) {
			cout << "Reading ptx file " << infile << " to memory..." << endl;

			string line;
			getline(ifs, line);
			const int width = stoi(line);
			getline(ifs, line);
			const int height = stoi(line);

			cout << "Image dimensions: " << width <<  " x " << height << endl;
		
			// skips header lines
			for(int i=0; i<8; i++) {
				getline(ifs, line);
			}

			unsigned int total_lines = width * height;
			vector<unsigned char> buf(4*total_lines);
			vector<double> xyz(3*total_lines);
			vector<double> depth(total_lines);

			double cx, cy, cz, intensity;
			int r, g, b;
			double dist, maxDist = -100, minDist = 100;

			double distsum = 0;

			double zmax = -100, zmin = 100;
		
			for (int x = 0; x < width; ++x) {
				cout << '\r' << x << '/'<< width;
				for (int y = height-1; y >= 0; --y) {

					ifs >> cx >> cy >> cz >> intensity >> r >> g >> b;
					dist = sqrt(cx*cx + cy*cy + cz*cz);
					distsum += dist;

					// need to compute index because ptx is ordered by columns and png is ordered by rows
					int index = y*width+x;

					buf[4*index] = static_cast<unsigned char>(r);
					buf[4*index+1] = static_cast<unsigned char>(g);
					buf[4*index+2] = static_cast<unsigned char>(b);
					buf[4*index+3] = 255;

					xyz[3*index] = cx;
					xyz[3*index+1] = cy;
					xyz[3*index+2] = cz;

					depth[index] = dist < 1e-6 ? 0 : dist;

					cz /= dist;

					zmax = max(cz, zmax);
					zmin = min(cz, zmin);
					maxDist = max(dist, maxDist);
					if (dist > 1e-6) {
						minDist = min(dist, minDist);
					}
				}
			}
		

			cout << endl;

			cout << "zmax: " << zmax << endl << "zmin: " << zmin << endl;
			cout << "maxdist: " << maxDist << endl << "mindist: " << minDist << endl;
			cout << "avgdist: " << (distsum/(width*height)) << endl;

			maxDist = 3;

			int response = 0;
			while(response != -1) {
				response = writeImage(width, height, buf, xyz, depth);
			}

			ifs.close();
		}
	} else if (infile.substr(infile.find_last_of('.')+1) == "png") {

		cout << "Reading png file " << infile << " to memory..." << endl;

		vector<unsigned char> png;
		vector<unsigned char> buf;
		unsigned width, height;

		lodepng::load_file(png, infile);
		unsigned error = lodepng::decode(buf, width, height, png);

		cout << "Image dimensions: " << width << "x" << height << endl;

		int response = 0;
		while(response != -1) {
			response = writeImage(width, height, buf);
		}
	}

	return 0;
}

// write an image file based on user input
int writeImage(int width, int height, vector<unsigned char>& buf, vector<double>& xyz, vector<double>& depth, double minDist, double maxDist) {
	string input = "";
	cout << "Enter params (or help): ";
	getline(cin, input);

	cout << endl;

	if (input == "help" || input == "-h" || input == "--help" || input == "") {
		cout << "Params: outfile width height theta phi verticalFOV\n  outfile: output file (*.png)\n  width: width of the output image\n  height: height of the output image"
				<< "\n  theta: theta angle in radians corresponding to image center\n  phi: phi angle in radians corresponding to image center\n  verticalFOV: vertical field of view in radians" << endl;
		return 0;
	}

	if (input == "quit" || input == "exit") {
		return -1;
	}

	stringstream ss(input);

	string outfile;
	int imageWidth, imageHeight;
	double theta, phi, fovy;

	ss >> outfile >> imageWidth >> imageHeight >> theta >> phi >> fovy;
			
	Vector3d w(cos(theta)*cos(phi), sin(theta)*cos(phi), sin(phi));
	Vector3d u = w.cross(Vector3d(0, 0, 1)).normalized();
	Vector3d v = w.cross(u).normalized();

	double vectorLength = 2 * tan(fovy / 2) / imageHeight;
	u *= vectorLength;
	v *= vectorLength;

	Vector3d origin = w - imageWidth/2 * u - imageHeight/2 * v;

	vector<unsigned char> outBuf(4*imageWidth*imageHeight);

	vector<double> xyzBuf;
	if (xyz.size() > 0) {
		xyzBuf.resize(3*imageWidth*imageHeight);
	}

	vector<unsigned char> depthBuf;
	if (depth.size() > 0) {
		depthBuf.resize(4*imageWidth*imageHeight);
	}

	int imageIndex = 0, xyzIndex = 0;
	for(int y=0; y<imageHeight; ++y) {
		for(int x=0; x<imageWidth; ++x) {

			Vector3d point = origin + x * u + y * v;
			point.normalize();

			double cx, cy;
			vector<double> blend;
			int index = getIndex(point, width, height, blend);

			double dx = cx - (int) cx;
			double dy = cy - (int) cy;

			outBuf[imageIndex] = static_cast<unsigned char>(blend[0] * buf[4*index] + blend[1] * buf[4*(index+1)] + blend[2] * buf[4*(index+width)] + blend[3] * buf[4*(index+width+1)]);
			outBuf[imageIndex+1] = static_cast<unsigned char>(blend[0] * buf[4*index+1] + blend[1] * buf[4*(index+1)+1] + blend[2] * buf[4*(index+width)+1] + blend[3] * buf[4*(index+width+1)+1]);
			outBuf[imageIndex+2] = static_cast<unsigned char>(blend[0] * buf[4*index+2] + blend[1] * buf[4*(index+1)+2] + blend[2] * buf[4*(index+width)+2] + blend[3] * buf[4*(index+width+1)+2]);
			outBuf[imageIndex+3] = static_cast<unsigned char>(blend[0] * buf[4*index+3] + blend[1] * buf[4*(index+1)+3] + blend[2] * buf[4*(index+width)+3] + blend[3] * buf[4*(index+width+1)+3]);

			if (xyz.size() > 0) {
				xyzBuf[xyzIndex] = xyz[3*index];
				xyzBuf[xyzIndex+1] = xyz[3*index+1];
				xyzBuf[xyzIndex+2] = xyz[3*index+2];
			}

			if (depth.size() > 0) {
				int depthColor = (int) ((1 - min((depth[index] - minDist) / (maxDist - minDist), 1.0)) * 255);
				depthBuf[imageIndex] = static_cast<unsigned char>(depthColor);
				depthBuf[imageIndex+1] = 0;
				depthBuf[imageIndex+2] = 0;
				depthBuf[imageIndex+3] = 255;
			}

			xyzIndex += 3;
			imageIndex += 4;

		}
	}

	vector<unsigned char> png;
	unsigned int error = lodepng::encode(png, outBuf, imageWidth, imageHeight);

	if (!error) {
		lodepng::save_file(png, outfile);
		cout << "SUCCESS: saved image file " << outfile << endl;
	} else {
		cout << "ERROR: error while encoding, this may be due to an incorrect specified number of pixels" << endl;
		cout << "Expected " << imageWidth*imageHeight << " pixels and found " << outBuf.size()/4 << endl;
	}

	if (xyz.size() > 0) {
		ofstream out(outfile.substr(0, outfile.length()-4) + ".obj");
		for(int i=0; i<imageWidth*imageHeight; i++) {
			out << "v " << xyzBuf[3*i] << " " << xyzBuf[3*i+1] << " " << xyzBuf[3*i+2] << " " 
				<< static_cast<double>(outBuf[4*i])/255.0 << " " << static_cast<double>(outBuf[4*i+1])/255.0 << " " << static_cast<double>(outBuf[4*i+2])/255.0 << endl;
		}

		cout << "SUCCESS: saved obj file " << outfile.substr(0, outfile.length()-4) + ".obj" << endl;

		out.close();
	}

	if (depth.size() > 0) {
		vector<unsigned char> depthpng;
		error = lodepng::encode(depthpng, depthBuf, imageWidth, imageHeight);

		if (!error) {
			lodepng::save_file(depthpng, "depth_" + outfile);
			cout << "SUCCESS: saved depth image file depth_" << outfile << endl;
		} else {
			cout << "ERROR: error while encoding depth image file, this may be due to an incorrect specified number of pixels" << endl;
			cout << "Expected " << imageWidth*imageHeight << " pixels and found " << depthBuf.size()/4 << endl;
		}
	}

	return 0;
}

// gets the index of the pixel based on a given vector
int getIndex(Vector3d& point, int width, int height, vector<double>& blend, bool debug) {

	double theta, phi;

	theta = atan(point[1] / point[0]);
	phi = asin(point[2] / point.norm());

	if (point[0] < 0) {
		theta += M_PI;
	} else if (point[1] < 0) {
		theta += 2*M_PI;
	}

	double x = (1 - (theta - THETA_MIN) / (THETA_MAX - THETA_MIN)) * width;
	double y = ((1 - (phi - PHI_MIN) / (PHI_MAX - PHI_MIN)) * height);

	double dx = x - (int) x;
	double dy = y - (int) y;

	// approximate
	blend.push_back((1-dx)*(1-dy));
	blend.push_back(dx * (1-dy));
	blend.push_back((1-dx) * dy);
	blend.push_back(dx * dy);

	int index = getIndex(width, height, x, y);

	if (debug) {
		cout << "XYZ: " << point[0] << " " << point[1] << " " << point[2] << endl;
		cout << "ANGLES: " << theta << " " << phi << " | " << index << endl;
	}

	return index;
}

int getIndex(int width, int height, double x, double y) {
	int index = ((int) x + (int) y * width) % (width * height); 

	if (index < 0) {
		index += width * height;
	}

	return index;
}