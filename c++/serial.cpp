#include "pythongrid_int.h"
#include <iostream>
#include <sstream>
#include <cstring>
int main()
{
double v=40;
std::ostringstream vname;
vname<<v;
double f [6]={0.0,0.2,0.4,0.6,0.8,1.0};
double accre [11]={0.0,0.6,1.2,1.8,2.4,3.0,3.6,4.2,4.8,5.4,6.0};
for (int fi=0; fi<6;fi++){
	for (int aci=0;aci<11;aci++){
                std::ostringstream fname;
                std::ostringstream acname;
                fname<<f[fi];
                acname<<accre[aci];
		std::string outname="../cubes/cserial/v"+vname.str()+"_f"+fname.str()+"_ac"+acname.str()+".fits";
                char* outnamechar=new char[outname.length() + 1];
                strcpy(outnamechar,outname.c_str());
                std::cout<<v<<" "<<f[fi]<<" "<<accre[aci]<<" "<<outnamechar<<std::endl;
                double re=grid((float)v,(float)f[fi],(float)accre[aci],outnamechar);

}
}










return 0;
}
