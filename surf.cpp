#include <queue> //in bfs
#include "surf.h"
using namespace std;

//constructor
SurfaceStruct::SurfaceStruct(const GridData &grid, double gamma):grid(grid){
	phi = matrix<double>(grid.nx[0],grid.nx[1],grid.nx[2]);	
	set_field<double>(grid,phi,1.0);//initialize as outside
}

//default constructor 
SurfaceStruct::SurfaceStruct(){}

//copy constructor
SurfaceStruct::SurfaceStruct(const SurfaceStruct& surf){
	grid = surf.grid;
	phi = matrix<double>(grid.nx[0],grid.nx[1],grid.nx[2]);	
	copy_field<double>(grid, surf.phi, phi );
}
//destructor
SurfaceStruct::~SurfaceStruct(){
	free_matrix<double>(phi, grid.nx[0], grid.nx[1], grid.nx[2]);
}

//swap
void swap(SurfaceStruct &first, SurfaceStruct &second){
	std::swap(first.p, second.p);
	std::swap(first.grid, second.grid);
	std::swap(first.phi, second.phi);
}
// assign 
SurfaceStruct& SurfaceStruct::operator = (SurfaceStruct rhs){
	swap(*this,rhs);
	return *this;
}

// read surface from file
void SurfaceStruct::read(string file){

	ifstream infile(file);
	if(infile.fail()){
		printf ("Error opening file");
        exit(1);
    }
	//read data from surffile
	for (int i=0;i<=grid.nx[0];i++){
		for (int j=0;j<=grid.nx[1];j++){
			for (int k=0;k<=grid.nx[2];k++){
					infile>>phi[i][j][k];
			}
		}
	}
	double temp;
	cout << "[remaining data]\n";
	while (infile >> temp) 
	{
		cout << temp << endl;
	}
	if (infile.eof())// check for EOF
		cout << "[EoF reached]\n";

	infile.close();
}


void SurfaceStruct::flip(int i, int j, int k){
	phi[i][j][k] = -phi[i][j][k];
}

void SurfaceStruct::flip(Sub sub){
	phi[sub[0]][sub[1]][sub[2]] = -phi[sub[0]][sub[1]][sub[2]];
}

void SurfaceStruct::flip(int ind){
	Sub sub = ind2sub(grid, ind);
	flip(sub);
}



bool SurfaceStruct::is_inside(Sub sub){
	return phi[sub[0]][sub[1]][sub[2]] < 0.0;
}

bool SurfaceStruct::is_outside(Sub sub){
	return phi[sub[0]][sub[1]][sub[2]] > 0.0;
}


// get direct neighbor subscript, up down left right
// might be out of bound
vector<Sub> SurfaceStruct::get_nbrs(Sub sub){
	vector<Sub> v;
	for(int i : {0,1,2}){
		for(int s : {-1,1}){
			Sub nbr = sub;
			nbr[i] = sub[i] + s;
			v.push_back(nbr);	
			
		}
	}
	return v;
}



// look at direct neighbor, check if different sign
// or neighbor is out of bound, self is inside
bool SurfaceStruct::near_interface(Sub sub){
	vector<Sub> nbrs = get_nbrs(sub);

	for( auto nbr : nbrs){
		if(outofbound(nbr, grid)){
			return phi[sub[0]][sub[1]][sub[2]]<0?true:false;
		} else{
			if (phi[nbr[0]][nbr[1]][nbr[2]] * phi[sub[0]][sub[1]][sub[2]] <= 0.0){
				// if change sign
				return true;	
			}	
		}
		
	}
	return false;
}

