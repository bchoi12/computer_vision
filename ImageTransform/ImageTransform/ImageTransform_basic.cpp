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

int getIndex(Vector3d& point, int width, int height, bool debug = false);

int main(int argc, char* argv[])
{
	if (argc < 1) {
		cout << "Usage: ImageTransform infile\n  infile: the input file (*.ptx)" << endl;
		exit(1);
	}

	string infile = argv[1];

	ifstream ifs;
	ifs.open(infile);

	if (ifs.is_open()) {
		cout << "Reading " << infile << " to memory..." << endl;

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

		double cx, cy, cz, intensity;
		int r, g, b;

		double zmax = -100, zmin = 100;
		
		for (int x = 0; x < width; ++x) {
			cout << '\r' << x << '/'<< width;
			for (int y = height-1; y >= 0; --y) {

				ifs >> cx >> cy >> cz >> intensity >> r >> g >> b;

				// need to compute index because ptx is ordered by columns and png is ordered by rows
				int index = 4*(y*width+x);

				buf[index] = static_cast<unsigned char>(r);
				buf[index+1] = static_cast<unsigned char>(g);
				buf[index+2] = static_cast<unsigned char>(b);
				buf[index+3] = 255;

				cz /= sqrt(cx*cx + cy*cy + cz*cz);

				if (cz > zmax) {
					zmax = cz;
				}
				if (cz < zmin) {
					zmin = cz;
				}
			}
		}
		

		cout << endl;

		cout << "zmax: " << zmax << endl << "zmin: " << zmin << endl;

		while(true) {
			string input = "";
			cout << "Enter params (or help): ";
			getline(cin, input);

			cout << endl;

			if (input == "help" || input == "-h" || input == "--help" || input == "") {
				cout << "Params: outfile width height theta phi verticalFOV\n  outfile: output file (*.png)\n  width: width of the output image\n  height: height of the output image"
					 << "\n  theta: theta angle in radians corresponding to image center\n  phi: phi angle in radians corresponding to image center\n  verticalFOV: vertical field of view in radians" << endl;
				continue;
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
			int imageIndex = 0;
			for(int y=0; y<imageHeight; ++y) {
				for(int x=0; x<imageWidth; ++x) {

					Vector3d point = origin + x * u + y * v;
					point.normalize();

					int index = getIndex(point, width, height);

					outBuf[imageIndex] = buf[index];
					outBuf[imageIndex+1] = buf[index+1];
					outBuf[imageIndex+2] = buf[index+2];
					outBuf[imageIndex+3] = buf[index+3];

					imageIndex += 4;

				}
			}

			vector<unsigned char> png;
			lodepng::State state;
			unsigned int error = lodepng::encode(png, outBuf, imageWidth, imageHeight, state);

			if (!error) {
				lodepng::save_file(png, outfile);
				cout << "SUCCESS: saved image file " << outfile << endl;
			} else {
				cout << "ERROR: error while encoding, this may be due to an incorrect specified number of pixels" << endl;
			}

		}
		ifs.close();
	}

	return 0;
}

// gets the index of the pixel based on a given vector
int getIndex(Vector3d& point, int width, int height, bool debug) {

	double theta, phi;

	theta = atan(point[1] / point[0]);
	phi = asin(point[2] / point.norm());

	if (point[0] < 0) {
		theta += M_PI;
	} else if (point[1] < 0) {
		theta += 2*M_PI;
	}

	int index = (int) ((1 - (theta - THETA_MIN) / (THETA_MAX - THETA_MIN)) * width + (int) ((1 - (phi - PHI_MIN) / (PHI_MAX - PHI_MIN)) * height) * width) % (width * height); 

	if (index < 0) {
		cout << "Warning: field of view may be too large" << endl;
		index += width * height;
	}

	if (debug) {
		cout << "XYZ: " << point[0] << " " << point[1] << " " << point[2] << endl;
		cout << "ANGLES: " << theta << " " << phi << " | " << index << endl;
	}

	return 4*index;
}

