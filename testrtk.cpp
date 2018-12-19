//
// Created by Haoliu on 8/1/18.
//

#include "GNSS.h"
#include <stdio.h>
#include "EphemSp3.h"
using namespace std;
using namespace Eigen;

void WriteSols(SolutionDeque sols,string saveName){
    string path = "../log/"+saveName+".txt";
    FILE *fp = fopen(path.data(),"w");
    for(Solution sl:sols){
//       fprintf(fp,"%.4f,%.4f,%.4f,%.4f\n",sl.time.tow,sl.xyz(0),sl.xyz(1),sl.xyz(2));
        fprintf(fp,"%.4f,%.8f,%.8f,%.8f\n",sl.time.tod,sl.lla(0)*R2D,sl.lla(1)*R2D,sl.lla(2));
    }
    fclose(fp);
}
    auto split = [](char* str,char c,vector<string> &result){
        char temp[16],*p = str;
        int k=0,n=0;
        while ((*p)&&14!=n){
            if(*p==c){
                temp[k]=0;
                k=0;n++;
                result.push_back(temp);
                memset(temp,0,16);
            } else{
                temp[k++]=*p;
            }
            p++;
        }
    };

int main(){
    FILE* fp = fopen("/home/alan/Downloads/novatel_inspvax.txt","r");
    FILE* fpw = fopen("../log/1219_07_41xyzNVT.txt","wb");
    char buffer[256];
    while (fgets(buffer, sizeof(buffer),fp)){
        vector<string> strs;
        split(buffer,',',strs);
        double time = atof(strs[0].c_str());
        printf("tiem= %f,%s\n", time,strs[0].c_str());
        fprintf(fpw,"%.4f",time/1000-86400*3);
        for (int i = 1; i < strs.size(); ++i) {
            fprintf(fpw,",%s",strs[i].c_str());
        }
        fprintf(fpw,"\n");
    }
    fclose(fp);fclose(fpw);

    double dd = 0.000001;
    Vector3d xyz_,xyzw,xyzj;
    Vector3d lla_(40*D2R,116.3*D2R,60);
    Vector3d llaw((40+dd)*D2R,116.3*D2R,60);
    Vector3d llaj(40.0*D2R,(116.3+dd)*D2R,60);

    LLA2XYZ(lla_,xyz_);
    LLA2XYZ(llaw,xyzw);
    LLA2XYZ(llaj,xyzj);
    ShowV3(xyz_,"xyz");
    ShowV3(xyzw,"xyzw");
    ShowV3(xyzj,"xyzj");
    printf("dn,dw %f,%f\n",(xyzw-xyz_).norm(),(xyzj-xyz_).norm() );
    return 0;
}
