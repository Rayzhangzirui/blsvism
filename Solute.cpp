#include "Solute.h"
#include "util.h"
#include <sstream>
#include <iomanip>
#include <algorithm>

//Effects of lengthscales and attractions on the collapse of hydrophobic polymers in water. Results and Discussion.
double c25sigma = 3.73;//LJss sigma; AA
double c25eps = 0.5856;//LJss eps: kj/mol

// read atom from file as matrix
vector<vector<double>> readatom(std::string filename){
	vector<vector<double>> mx;
	ifstream indata(filename,ios::in);
	
	if(!indata){
		cerr<<"open file failed";
		exit(1);
	}

	string line;
	int k = 0; //atom counter
	while(getline(indata,line)){
		if(line[0] =='#' || line.empty()){
			continue;
		}
		std::replace(line.begin(), line.end(), ',' , ' ');
		istringstream iss(line);
		vector<double> row;
		double tmp;
		while(iss >> tmp){
			row.push_back(tmp);
		}
		
		mx.push_back(row);

	}
	return mx;
}

SoluteStruct CreateC25(int n){
	double pos[25][3]=
    {
        {-1.11576,    -0.0676292,    -0.840911},
        {-0.0252184,    0.848302,    -0.283428},
        {-0.00760816,    2.18755,    -1.01787},
        {-0.58358,    3.18215,    -0.067917},
        {-0.657728,    4.41703,    -0.916446},
        {-2.05408,    5.00679,    -0.714063},
        {-2.969,    4.50911,    -1.77389},
        {-3.97513,    3.55555,    -1.13188},
        {-3.68004,    2.11315,    -1.40136},
        {-4.77188,    1.35306,    -0.679467},
        {-4.08921,    0.379089,    0.279966},
        {-4.00142,    0.899636,    1.74868},
        {-2.8411,    1.85559,    1.89361},
        {-3.34449,    2.85297,    2.86604},
        {-3.31371,    4.21235,    2.16758},
        {-2.49824,    5.16136,    3.02558},
        {-1.1412,    5.59752,    2.54612},
        {-0.0691338,    5.06892,    3.50928},
        {0.649921,    3.95055,    2.78243},
        {0.384614,    2.72983,    3.64455},
        {0.453911,    1.43355,    2.86057},
        {-0.451216,    0.325864,    3.44676},
        {-0.350434,    -0.880567,    2.60273},
        {-1.72801,    -1.51623,    2.50829},
        {-1.80904,    -2.2657,    1.20567},
    };

    SoluteStruct sol(n);

    for(int i = 0; i < n; i++){
    	for(int j = 0; j < 3; j++){
    		sol.d[i][j] = pos[i][j];
    	}
    	sol.sigma[i] = c25sigma;
    	sol.epsilon[i] = c25eps;
		sol.Q[i] = 0;
    }
    return sol;
}

SoluteStruct SolFromFile(string filename){
	if (filename.empty()){
		cout<<"empty filename, use default c25" <<endl;
		return CreateC25(25);
	} else if (filename=="methane") {
		return CreateMethane();
	}

	vector<vector<double>> mx = readatom(filename);
	int  n = mx.size();
	int  m = mx[0].size();
	cout<<"row " << n << ", column " << m<<endl;

	SoluteStruct sol(n);

	if ( m == 3) {
		cout << "input file coordinate only" <<endl;
		for(int i = 0; i < n; i++){
	    	for(int j = 0; j < 3; j++){
	    		sol.d[i][j] = mx[i][j];
	    	}
	    	sol.sigma[i] = c25sigma;
	    	sol.epsilon[i] = c25eps;
			sol.Q[i] = 0;
	    }		
	} else if ( m == 7) {
		cout << "input file all" <<endl;
		for(int i = 0; i < n; i++){
	    	for(int j = 1; j <= 3; j++){
	    		sol.d[i][j-1] = mx[i][j];
	    	}
	    	
			sol.sigma[i] = mx[i][4];
			sol.epsilon[i] = mx[i][5];
			sol.Q[i] = mx[i][6];
	    }		
	} else {
		cerr << "unknow matrix size"<<endl;
		exit(1);
	}
	sol.Print(3);
	return sol;
}



// make solute consists of n atom at (0,0,0), (1,0,0), ...
SoluteStruct CreateMethane(int n, double charge, double eps, double sigma){
	SoluteStruct sol(n);
	for (int i = 0; i < n; i ++){
		sol.d[i][0] = 0 + i;
		sol.d[i][1] = 0;
		sol.d[i][2] = 0;
		sol.epsilon[i] = eps;
		sol.sigma[i] = sigma;
		sol.Q[i] = charge;	
	}
	
	return sol;
}

// empty constructor
SoluteStruct::SoluteStruct(){}
//constructor intialize grid
SoluteStruct::SoluteStruct(int n){
	dnum = n;
	d = matrix<double>(dnum-1,2); //matrix of size dnum x 3
	index = new int[dnum];
	epsilon = new double[dnum];
	sigma = new double[dnum];
	Q = new double[dnum];
}


//copy constructor, delegate to constructor
SoluteStruct::SoluteStruct(const SoluteStruct& s):SoluteStruct(s.dnum){
	//copy data
	copy(s.index, s.index+dnum, index);
	copy(s.epsilon, s.epsilon+dnum, epsilon);
	copy(s.sigma, s.sigma+dnum, sigma);
	copy(s.Q, s.Q+dnum, Q);

	for (int i = 0; i<s.dnum; i++){
		copy(s.d[i],s.d[i]+3, d[i]);
	}
}

//destructor
SoluteStruct::~SoluteStruct(){
	free_matrix<double>(d, dnum-1, 3-1);
	delete [] index;
	delete [] epsilon;
	delete [] sigma;
	delete [] Q;

}

// get centroid of solute
va::Vector SoluteStruct::centroid() const{
	double temp_x = 0;
	double temp_y = 0;
	double temp_z = 0;
	for(int i = 0; i < dnum; i++){
		temp_x += d[i][0];
		temp_y += d[i][1];
		temp_z += d[i][2];
	}
	temp_x = temp_x / dnum;
	temp_y = temp_y / dnum;
	temp_z = temp_z / dnum;
	return va::Vector(temp_x,temp_y,temp_z);
}


double SoluteStruct::RadiusGyration() const{
	va::Vector c = centroid();

	double Rg=0.;

	for(int i=0; i<dnum; i++)
	{
		Rg+=(d[i][0]-c.x())*(d[i][0]-c.x())+
			(d[i][1]-c.y())*(d[i][1]-c.y())+
			(d[i][2]-c.z())*(d[i][2]-c.z());
	}
	return sqrt(Rg/dnum);
}


void SoluteStruct::Print(int n) const{
	for (int i = 0; i < min(n,dnum); i++){
		cout<< i << " ("<< setprecision(3)<< left << d[i][0] << ", "<<d[i][1] << ", "<<d[i][2]<<") "
			<< sigma[i] <<", " << epsilon[i] << ", " << Q[i] <<endl;
	}

}


void SoluteStruct::Write(string fname) const{
	ofstream fo(fname);
	for (int i = 0; i < dnum; i++){
		fo<< fixed << setprecision(3)<< d[i][0] << " "<<d[i][1] << " "<<d[i][2]<<endl;
	}

}
// rotate solute
void SoluteStruct::rotate(va::Rotation R){
		#ifdef DEBUG
		cout<<"SoluteStruct::rotate "<<R<<endl;
		#endif
		va::rmatrix q_matrix = to_rmatrix( R );

		va::Vector c = centroid();
		double center_x = c.x();
		double center_y = c.y();
		double center_z = c.z();
		#pragma parallel for
		for(int i = 0; i < dnum; i++){
			double x =d[i][0] - center_x;
			double y =d[i][1] - center_y;
			double z =d[i][2] - center_z;

			d[i][0] = q_matrix.a11 * x + q_matrix.a12 * y + q_matrix.a13 * z + center_x;
			d[i][1] = q_matrix.a21 * x + q_matrix.a22 * y + q_matrix.a23 * z + center_y;
			d[i][2] = q_matrix.a31 * x + q_matrix.a32 * y + q_matrix.a33 * z + center_z;
		}
}

void SoluteStruct::rotate(va::Vector u, double th){
		va::Rotation R( u, th );
		rotate(R);
}

void SoluteStruct::rotate(double ax, double ay, double az, double th){
	rotate(va::Vector(ax,ay,az),th);
}

// translate solute
void SoluteStruct::translate(va::Vector v){
	#pragma parallel for
	for (int i=0; i<dnum; i++){
		d[i][0] = d[i][0] + v.x();
		d[i][1] = d[i][1] + v.y();
		d[i][2] = d[i][2] + v.z();
	}
}


void swap(SoluteStruct &first, SoluteStruct &second){
	
	std::swap(first.dnum, second.dnum);
	std::swap(first.index, second.index);
	std::swap(first.d, second.d);
	std::swap(first.epsilon, second.epsilon);
	std::swap(first.sigma, second.sigma);
	std::swap(first.rho0, second.rho0);
	std::swap(first.estat, second.estat);
	std::swap(first.epsilonex, second.epsilonex);
	std::swap(first.epsilonin, second.epsilonin);
	std::swap(first.epsilon0, second.epsilon0);	
};

// assignment =
SoluteStruct& SoluteStruct::operator = (SoluteStruct other){
	swap(*this, other);
	return *this;
};

//c2c distance
double SoluteStruct::c2c_dist(const SoluteStruct &p2) const{
	va::Vector v = p2.centroid() - this->centroid();
	return v.norm();
}

//minimimum distance
double SoluteStruct::min_pairwise_dist(const SoluteStruct &p2) const {
	double rmin = 1e12;
	double r;
	for (int i=0; i<dnum; i++){
		for (int j=0; j<p2.dnum; j++){
			r = sqrt(pow(d[i][0] - p2.d[i][1], 2)+
					pow(d[i][1] - p2.d[i][1], 2)+
					pow(d[i][2] - p2.d[i][1], 2));
			rmin = (r<rmin)?r:rmin;
		}
	}
	return rmin;
};

// https://stackoverflow.com/questions/1657883/variable-number-of-arguments-in-c
// get center of 2 solute
va::Vector SoluteStruct::get_center(std::initializer_list<SoluteStruct> solutes){
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;

	for(auto &solute : solutes)
	{
		for(int i = 0; i < solute.dnum; i++){
			x.push_back(solute.d[i][0]);
			y.push_back(solute.d[i][1]);
			z.push_back(solute.d[i][2]);	
		}
		
	}
	double min_x = *std::min_element(x.begin(), x.end());
	double min_y = *std::min_element(y.begin(), y.end());
	double min_z = *std::min_element(z.begin(), z.end());

	double max_x = *std::max_element(x.begin(), x.end());
	double max_y = *std::max_element(y.begin(), y.end());
	double max_z = *std::max_element(z.begin(), z.end());

	va::Vector v((min_x + max_x)/2,(min_y + max_y)/2,(min_z + max_z)/2);

	cout<<"bound: " <<min_x << " " << max_x << " " << min_y << " " << max_y<< " " << min_z<< " " << max_z<<endl;
	return v;
}



// read data from file, start index, end index
void SoluteStruct::readdata(std::string filename, int start, int end){
	if(end-start+1!=dnum){
		cerr<<"atom numbers does not match. dnum = "<<dnum<<", file = "<<end-start+1<<endl;
		exit(1);
	}
	ifstream indata(filename,ios::in);
	if(!indata){
		cerr<<"open file failed";
		exit(1);
	}
	string line;
	int k = 0; //atom counter
	int ind;
	double x,y,z,s,e,q;
	while(getline(indata,line)){
		if(line[0]!='#'){
			//skip comment line with #
			istringstream iss(line);
			iss>>ind>>x>>y>>z>>s>>e>>q;
			if(start <=ind && ind <= end){
				index[k] = ind;
				d[k][0] = x;
				d[k][1] = y;
				d[k][2] = z;
				sigma[k] = s;
				epsilon[k] = e;
				Q[k] = q;
				k++;
			}
		}
	}
}


// read data from file, start index, end index
void SoluteStruct::readpos(std::string filename){
	cout<<"read pos"<<endl;
	ifstream indata(filename,ios::in);
	if(!indata){
		cerr<<"open file failed";
		exit(1);
	}
	
	for (int i = 0; i < dnum ; i ++ ) {
		for (int j = 0; j < 3 ; j ++ ) {
			if (!(indata>>d[i][j])) {
				exit(1);
			}
		}
	}
	// cout<<"unread pos"<<endl;
	// cout << indata.rdbuf();	
}


// maximum coordinate in 
double SoluteStruct::MaxCoord() const{
	double m = 0.0;

	for (int i=0; i<dnum; i++){
		for (int j=0; j<3; j++){
			if( abs(d[i][j]) > m ){
				m = abs(d[i][j]);
			}
		}
	}
	return m;
}


void SoluteStruct::Centering(){
	va::Vector center = get_center({*this});
    translate(center*(-1.0));
}