//
//
//#include <complex.h>
//#include <typeinfo>
//#include <iostream>
//#include <stdlib.h>
//#include <string.h>
//
//const int FLT = 0 ;
//const int DBL=1;
//const int CPLXFLT =2 ;
//const int CPLXDBL= 3;
//
//template <typename _Tp>
//int gettype(){
//
//        double dbl = 0.;
//        float flt = 0.;
//        complex double cdbl = 0. +I;
//        complex float cflt = 0. + I ;
//
//
//        if(typeid(_Tp) == typeid(flt))
//                return FLT;
//        if(typeid(_Tp) == typeid(dbl))
//                return DBL;
//        if(typeid(_Tp) == typeid(cflt))
//                return CPLXFLT;
//        if(typeid(_Tp) == typeid(cdbl))
//                return CPLXDBL;
//
//        return -1;
//
//}
//
//
//
//
//template<typename _Tp>
//void* bytewisecopy(_Tp* val,int type, int length)
//{
//        unsigned int bytesize = 0 ;
//
//        if (type == FLT)
//                bytesize = sizeof(float);
//        else if (type == DBL)
//                bytesize = sizeof(double);
//        else if (type == CPLXFLT)
//                bytesize = sizeof(complex float);
//        else if (type == CPLXDBL)
//                bytesize = sizeof(complex double);
//	std::cout<<" bytesize  = " << bytesize << "\n";
//
//        void* result = NULL;
//        if(bytesize > 0)
//        {
//                //copy only the relevant bytes!!!
//                result = (void*)malloc(length*bytesize);
//        	void * vval = (void*) val;
//		memcpy(result,vval,bytesize*length);
//
//
// 		std::cout<< "\n";
//		std::cout<<"val = " <<*val << "\n";
//		std::cout<<"res = " << *(_Tp*)(result)<<"\n";
//	}
//
//        return result;
//}
//
//
//
//int main(int argc, char** argv)
//{
//
//	double d = 1;
//	float f =  1.;
//	complex double cd = 1. + I;
//	complex float cf = 1. + I;
//
//	void* c = bytewisecopy<double>(&d,DBL,1);
//
//	float cval = *((float*)&cf);
//	std::cout << "Double to FLoat -> " << cval << "\n";
//	free(c);
//
//
//
//
//
//
//
//
//}
//
//
//
//
