#include "Mesh.hpp"

Mesh3D::Mesh3D(Vector3 domain, double oxgD, double egfD, double oBorder, double hCoarse, double hRefined, double deltaT){
    this->hCoarse = hCoarse;
    this->hRefined = hRefined;
    this->MatrixRefinedSize = 0;
    this->unityCoarse = Vector3i((domain/hCoarse)+1);
    this->unityRefined = Vector3i((domain/hRefined)+1);
    this->matrixSize = this->unityCoarse.x * this->unityCoarse.y * this->unityCoarse.z;
    this->uO = std::vector<double>(this->matrixSize, 0.0);
    this->nzNum = 8*4 + 4*(this->unityCoarse.x-2)*5 + 4*(this->unityCoarse.y-2)*5 + 4*(this->unityCoarse.z-2)*5 + 2*(this->unityCoarse.y-2)*(this->unityCoarse.z-2)*6 + 2*(this->unityCoarse.x-2)*(this->unityCoarse.z-2)*6 + 2*(this->unityCoarse.x-2)*(this->unityCoarse.y-2)*6 + (this->unityCoarse.x-2)*(this->unityCoarse.y-2)*(this->unityCoarse.z-2)*7;
    this->hReason = hCoarse/hRefined;
    this->deltaT = deltaT;

    this->sigma = oxgD*deltaT/pow(this->hCoarse, 2);

    this->B = new double[this->matrixSize];
    this->rowPtr = std::vector<int>((this->matrixSize)+1, 0.0);
    this->colInd = new int[this->nzNum];
    this->val = new double[this->nzNum];

    this->pos = std::vector<Vector3>(this->matrixSize, Vector3());

    for (size_t i = 0; i < std::max(this->unityRefined.x, std::max(this->unityRefined.y, this->unityRefined.z)); i++) {
        this->refined.push_back(Vector3(i*this->hRefined, i*this->hRefined, i*this->hRefined));
    }

    zera_crs(this->matrixSize, this->nzNum, &this->rowPtr[0], this->colInd, this->val);

    for (int i = 0; i < this->matrixSize; i++){
        this->pos[i].z = (i/(this->unityCoarse.z*this->unityCoarse.y))*this->hCoarse;
        this->pos[i].y = ((i%(this->unityCoarse.z*this->unityCoarse.y))/this->unityCoarse.z)*this->hCoarse;
        this->pos[i].x = ((i%(this->unityCoarse.z*this->unityCoarse.y))%this->unityCoarse.z)*this->hCoarse;
        for (int j = 0; j < this->matrixSize; j++){
            if (i==j){
                create_crs_ordened(i, j, &this->rowPtr[0], this->colInd, this->val, 0);
                continue;
            }
            if ( (this->pos[i].x!=0 && i-j==1) || (this->pos[i].x!=domain.x && j-i==1) ){
                create_crs_ordened(i, j, &this->rowPtr[0], this->colInd, this->val, -this->sigma);
                continue;
            }
            if ( (this->pos[i].y!=0 && i-j==this->unityCoarse.z) || (this->pos[i].y!=domain.y && j-i==this->unityCoarse.z) ){
                create_crs_ordened(i,j, &this->rowPtr[0], this->colInd, this->val, -this->sigma);
                continue;
            }
            if ( (i-j==(this->unityCoarse.z*this->unityCoarse.y)) || (j-i==(this->unityCoarse.z*this->unityCoarse.y)) ){
                create_crs_ordened(i,j, &this->rowPtr[0], this->colInd, this->val, -this->sigma);
            }
        }
    }

    //Matrix refined size of circulate domain
    Vector3 R;
	for (int i = 0; i < this->unityRefined.x; i++){
		for (int j = 0; j < this->unityRefined.y; j++){
			for (int k = 0; k < this->unityRefined.z; k++){
				R = Vector3 (pos[i].x-domain.x/2.0,pos[j].y-domain.y/2.0,pos[k].z-domain.z/2.0);
				if (sqrt(pow(R.x,2)+pow(R.y,2)+pow(R.z,2)) <= domain.x/2.0)
					this->MatrixRefinedSize++;
			}
		}
	}

    //Condiçãthis->unityCoarse.z de Contorno oxigênio
    for(int i=0;i<this->matrixSize;i++){
        this->B[i] = 0;
        //Dirichlet
        if (this->pos[i].x == 0 || this->pos[i].x == domain.x || this->pos[i].y == 0 || this->pos[i].y==domain.y || this->pos[i].z == 0 || this->pos[i].z==domain.z){
            this->B[i]= oBorder;
            for(int j = this->rowPtr[i]; j < this->rowPtr[i+1]; j++){
                if (this->colInd[j] == i){
                    this->val[j]=1.0;
                }
                else {
                    this->val[j]=0.0;
                }
            }
        }
    }
}


Mesh2D::Mesh2D(Vector3 domain, double oxgD, double egfD, double oBorder, double hCoarse, double hRefined, double deltaT){
    this->hCoarse = hCoarse;
    this->hRefined = hRefined;
    this->MatrixRefinedSize = 0;
    this->unityCoarse = Vector3i((domain.x/hCoarse)+1,(domain.y/hCoarse)+1,0.0);
    this->unityRefined = Vector3i((domain.x/hRefined)+1,(domain.y/hRefined)+1,0.0);
    this->matrixSize = this->unityCoarse.x * this->unityCoarse.y;
    this->uO = std::vector<double>(this->matrixSize, 0.0);

    this->nzNum = 4*3 + 2*(this->unityCoarse.x-2)*4 + 2*(this->unityCoarse.y-2)*4 + (this->unityCoarse.x-2)*(this->unityCoarse.y-2)*5;
    this->hReason = hCoarse/hRefined;
    this->deltaT = deltaT;

    this->sigma = deltaT*oxgD/pow(this->hCoarse, 2);

    this->B = new double[this->matrixSize];
    this->rowPtr = std::vector<int>((this->matrixSize)+1, 0.0);
    this->colInd = new int[this->nzNum];
    this->val = new double[this->nzNum];

    this->pos = std::vector<Vector3>(this->matrixSize, Vector3());

    for (size_t i = 0; i < std::max(this->unityRefined.x, this->unityRefined.y); i++) {
        this->refined.push_back(Vector3(i*this->hRefined, i*this->hRefined, 0.0));
    }

    zera_crs(this->matrixSize, this->nzNum, &this->rowPtr[0], this->colInd, this->val);

    for (int i = 0; i < this->matrixSize; i++){
        this->pos[i].y = (i/this->unityCoarse.x)*this->hCoarse;
        this->pos[i].x = (i%this->unityCoarse.x)*this->hCoarse;
        for (int j = 0; j < this->matrixSize; j++){
            if (i==j){
                create_crs_ordened(i, j, &this->rowPtr[0], this->colInd, this->val, 0);;
                continue;
            }
            if ( (this->pos[i].x!=0 && i-j==1) || (this->pos[i].x!=domain.x && j-i==1) ){
                create_crs_ordened(i, j, &this->rowPtr[0], this->colInd, this->val, -this->sigma);
                continue;
            }
            if ( (this->pos[i].y!=0 && i-j==this->unityCoarse.x) || (this->pos[i].y!=domain.y && j-i==this->unityCoarse.x) ){
                create_crs_ordened(i,j, &this->rowPtr[0], this->colInd, this->val, -this->sigma);
                continue;
            }
        }
    }

    //Matrix refined size of circulate domain
    Vector3 R;
	for (int i = 0; i < this->unityRefined.x; i++){
		for (int j = 0; j < this->unityRefined.y; j++){
			R = Vector3 (pos[i].x-domain.x/2.0,pos[j].y-domain.y/2.0,0.0);
			if (sqrt(pow(R.x,2)+pow(R.y,2)+pow(R.z,2)) <= domain.x/2.0)
				this->MatrixRefinedSize++;
		}
	}

    //Condiçãthis->unityCoarse.z de Contorno oxigênio
    for(int i=0;i<this->matrixSize;i++){
        this->B[i] = 0;
        //Neuman
        if (this->pos[i].x == 0 || this->pos[i].x == domain.x || this->pos[i].y == 0 || this->pos[i].y==domain.y ){
            this->B[i]= 0.0;
            for(int j = this->rowPtr[i]; j < this->rowPtr[i+1]; j++){
                if (this->colInd[j] == i){
                    this->val[j]=1.0;
                }
                else {
                    this->val[j]=0.0;
                }
                if (this->pos[i].x == 0 && this->colInd[j]-i == 1 && this->pos[i].y != 0 && this->pos[i].y != domain.y) this->val[j]=-1.0;
                if (this->pos[i].x == domain.x && i-this->colInd[j] == 1 && this->pos[i].y != 0 && this->pos[i].y != domain.y) this->val[j]=-1.0;
                if (this->pos[i].y == 0 && this->colInd[j]-i == this->unityCoarse.x) this->val[j]=-1.0;
                if (this->pos[i].y==domain.y && i-this->colInd[j] == this->unityCoarse.x) this->val[j]=-1.0;
            }
        }
    }
}
