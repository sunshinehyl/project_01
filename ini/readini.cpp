#include <io.h>
#include <direct.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include "INIReader.h"

using std::string;

void copyfile(const string& read_flie,string write_file){
    string temp;
    std::ifstream infile;
    std::ofstream outfile;

    infile.open(read_flie,std::ios::in);
    if(!infile)
        std::cout<<"The read-file open error"<<"\n";
    outfile.open(write_file,std::ios::out);
    if(!outfile)
        std::cout<<"The write-file open error"<<"\n";

    while (!infile.eof()){
        getline(infile, temp, '\n');
        outfile <<temp<<"\n";
    }
    outfile.close();
    infile.close();
}

void Split_String(const std::string& s, std::vector<std::string>& v)
{
  string temp="";
  for(int i;i<s.length();i++){
      if((s[i]==' ')&&(temp.length()>0)){
          v.push_back(temp);
          temp = "";
      }
      else if(s[i]!=' '){
          temp += s[i];
      }
  }
  v.push_back(temp);

}
//以特定的字符分割字符串
void SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
{
  std::string::size_type pos1, pos2;
  pos2 = s.find(c);
  pos1 = 0;
  while(std::string::npos != pos2)
  {
    v.push_back(s.substr(pos1, pos2-pos1));
 
    pos1 = pos2 + c.size();
    pos2 = s.find(c, pos1);
  }
  if(pos1 != s.length())
    v.push_back(s.substr(pos1));
}
//生成\CRWG_CSWG_VSIE_Tet.DAT文件
void Ctet_file(INIReader reader,const string& ctet_file){
    std::ofstream odat_CTet;
    odat_CTet.open(ctet_file,std::ios::out);
    if(!odat_CTet)
        std::cout<<"The write-file open error"<<"\n";
    odat_CTet<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<reader.GetReal("BACKGROUND_MEDIA", "EpsilonRe", 0)<<"    "
             <<reader.GetReal("BACKGROUND_MEDIA", "EpsilonIm", 0)<<"\t\t"<<'#'<<" BACKGROUND_MEDIA: EpsilonRe and EpsilonIm"<<"\n";
    odat_CTet<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<reader.GetReal("BACKGROUND_MEDIA", "MuRe", 0)<<"    "
             <<reader.GetReal("BACKGROUND_MEDIA", "MuIm", 0)<<"\t\t"<<'#'<<" BACKGROUND_MEDIA: MuRe and MuIm"<<"\n";

    odat_CTet.close();
}
//生成INCF_BiRCS.DAT文件
void Incf_file(INIReader reader,const string& incf_flie){
    std::vector<std::string> v_theta;
    std::vector<std::string> v_phi;
    string split_=":";
    std::ofstream odat_INCF;

    float a[3];
    float b[3];

    odat_INCF.open(incf_flie,std::ios::out);
    if(!odat_INCF)
        std::cout<<"The write-file open error"<<"\n";
    SplitString(reader.GetString("INCIDENT_WAVE", "Theta", "0:0:0"),v_theta,split_);
    
    for(std::vector<string>::size_type i = 0; i != v_theta.size(); ++i)
        a[i]=atof(v_theta[i].c_str());
    odat_INCF<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<"    "<<a[1]<<"\t\t\t"<<'#'
             <<"  Electric field amplitude of incident plane electromagnetic wave"<<"\n";
    odat_INCF<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<"    "<<a[2]<<"    "
             <<a[0]<<"\t\t"<<'#'<<"  Spherical coordinates of polarization direction /° TheE0,PhiE0"<<"\n";
    
    SplitString(reader.GetString("INCIDENT_WAVE", "Phi", "0:0:0"),v_phi,split_);
    for(std::vector<string>::size_type i = 0; i != v_phi.size(); ++i)
        b[i]=atof(v_phi[i].c_str());
    odat_INCF<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<"    "<<b[2]<<"    "
             <<b[0]<<"\t\t"<<'#'<<"  Spherical coordinates of propagation direction /°"<<"\n";
    odat_INCF.close();
}
//生成WFB_BiRCS.DAT文件
void Wfb_file(INIReader reader,const string& wfb_flie){
    std::vector<std::string> v_wfb;
    string split_=":";
    std::ofstream odat_Wfb;
    odat_Wfb.open(wfb_flie,std::ios::out);
    if(!odat_Wfb)
        std::cout<<"The write-file open error"<<"\n";
    //FREQUENCY
    int freq_num=0;
    SplitString(reader.GetString("FREQUENCY", "Value", "0:0:0"),v_wfb,split_);
    freq_num = (atof(v_wfb[2].c_str())-atof(v_wfb[0].c_str()))/atof(v_wfb[1].c_str());

    odat_Wfb<<0.5<<"\t\t"<<'#'<<"  Integral weight of CFIE [1.0]=EFIE;  [0.0]=MFIE\n";
    odat_Wfb<<atof(v_wfb[0].c_str())<<"\t\t"<<'#'<<" Freq0 Initial frequency of incident wave\n";
    odat_Wfb<<atof(v_wfb[1].c_str())<<"\t\t"<<'#'<<" Dfreq Calculate the interval between frequencies\n";
    odat_Wfb<<freq_num<<"\t\t"<<'#'<<" Numfreq Calculate the numeber of frequencies\n";

    //SLOVER-Preconditioner
    int Pre_num=3;
    string Pre_;
    Pre_ = reader.GetString("SOLVER", "Preconditioner", "no");
    if(Pre_=="diagonal"){
        Pre_num =1;
    }
    else if(Pre_=="normalized BF"){
        Pre_num =2;
    }
    else{
        Pre_num=3;
    }
    odat_Wfb<<Pre_num<<"    	# preconditioner：1->diagonal, 2->normalized BF, 3->no     #normalized BF produces very good results\n";
    //SLOVER-type
    int type_num=5;
    string type_;
    type_ = reader.GetString("SOLVER", "Type", "no");
    if(type_=="Gauss"){
        type_num =1;
    }
    else if(type_=="CGM"){
        type_num =2;
    }
    else if(type_=="BCGM"){
        type_num =3;
    }
    else if(type_=="GMRES"){
        type_num =4;
    }
    else if(type_=="LU"){
        type_num =5;
    }
    odat_Wfb<<type_num<<"    	# 1->Gauss, 2->CGM, 3->BCGM, 4->GMRES, 5->LU;\
             (notice:Gauss method is not suitable for calculation of RCS of single station / CGB or GMRES is not preferred for radiation problems)\n";
    //SLOVER-residual/Max_Iter
    odat_Wfb<<reader.GetReal("SOLVER", "Residual", 0)<<"	   		# iterative residual 0.01/0.001\n";
    odat_Wfb<<reader.GetReal("SOLVER", "Max_Iter", 0)<<"     			# max iterative number\n";
    odat_Wfb<<"40                 	# restart number of GMRES 40/30\n";
    //
    odat_Wfb<<"\n";
    //result-current
    int re_cu=0;
    string cu_;
    cu_ = reader.GetString("RESULT", "Current", "false");
    transform(cu_.begin(), cu_.end(), cu_.begin(), ::tolower);
    //判断current是true/false
    if(cu_=="true"){
        re_cu =0;
    }
    else if(cu_=="false"){
        re_cu =1;
    }
    odat_Wfb<<re_cu<<"       1     # Output_Ctype Output surface body (electric field) fluid in metal medium <1>y；<2>n\
              # !Voum_curr_type According to the field points in different coordinate systems， output electric field<1>Spherical coordinate；<2>Rectangular coordinate\n";
    odat_Wfb<<reader.GetReal("NEARFIELD", "NearFieldNumber", 0)<<"     100000             # Near_field, DL   !  1=>Yes；else=>No；else=>NEARFIELD/FARFIELD, DL Delineation area of NEARFIELD\n";
    odat_Wfb<<"2               # Total_Field !Under the excitation of point power：<1>Calculation： incident + scattering fields；<2>scattering\n";
    //频点
    if(freq_num>1)
        odat_Wfb<<"2              # Farfield_BiRCS  !  1=>BiRCS；2=>MonRCS；\n";
    else
        odat_Wfb<<"1              # Farfield_BiRCS  !  1=>BiRCS；2=>MonRCS；\n";
    
    odat_Wfb.close();
}
//生成elem.geo文件
int mid_elem(std::ifstream& infile, const string& out_file, int& tetrahedron_num, int& triangle_num){
    std::ofstream outfile_midelem ;

    outfile_midelem.open(out_file,std::ios::out);

    if((!outfile_midelem)||(!infile))
        return -6;

    string temp; //记录每一行的数据
    std::vector<std::string> v; 

    while(temp!="$Elements"){
        getline(infile, temp, '\n');
    }
    getline(infile, temp, '\n');
    getline(infile, temp, '\n');

    //将内容写到中间文件中
    while(temp!="$EndElements"){
        Split_String(temp,v);
        if(v.size()<2)
            std::cout<<"elements error";
        //填充四面体的数据
        if(v[1]=="4"){
            outfile_midelem<<"\t"<<v[v.size()-1];
            outfile_midelem<<"\t"<<v[0];
            for(std::vector<string>::size_type i = 2; i < 11; ++i)
                outfile_midelem<<"\t"<<v[i];
            outfile_midelem<<"\n\t"<<v[v.size()-2]<<"\n";
            tetrahedron_num++;   //计算tetrahedron_num和triangle_num的数量
        }
        //填充三角形的数据
        if(v[1]=="2"){
            triangle_num++;
            outfile_midelem<<"\t"<<2;
            outfile_midelem<<"\t"<<triangle_num;
            for(std::vector<string>::size_type i = 2; i < 8; ++i)
                outfile_midelem<<"\t"<<v[i];
            outfile_midelem<<"\t0"<<"\t0"<<"\t0"<<"\n"<<"\t0"<<"\n";
        }
        getline(infile, temp, '\n');
        v.clear();  
    }
    outfile_midelem.close();

    return 0;

}
int elem_geo(std::ifstream& infile, const string& out_file){
    string out_elem = out_file + "\\elem.geo";
    string out_midelem = out_file + "\\out_elem.geo";

    std::ofstream outfile_elem ;
    std::ifstream mid_file;
    
    string temp;

    int tetrahedron_num=0;
    int triangle_num = 0;
    
    mid_elem(infile,out_midelem,tetrahedron_num,triangle_num);
    
    mid_file.open(out_midelem.c_str(),std::ios::in);
    outfile_elem.open(out_elem.c_str(),std::ios::out);

    if((!outfile_elem)||(!mid_file))
        return -6;

    //加入四面体和三角形的数量
    outfile_elem <<"\t"<<tetrahedron_num<<"\n";
    outfile_elem <<"\t"<<triangle_num<<"\n";
    while (!mid_file.eof()){
        getline(mid_file, temp, '\n');
        outfile_elem <<temp<<"\n";
    }

    mid_file.close();
    outfile_elem.close();
    remove(out_midelem.c_str());

    return 0;
}
//生成elem.geo文件
void node_geo(std::ifstream& infile, std::vector<std::string>& geo_v, const string& out_file){
    string out_filename = out_file+ "\\node.geo";
    std::ofstream outfile_node ;

    outfile_node.open(out_filename.c_str(),std::ios::out);
    if((!outfile_node||(!infile)||(geo_v.size()<4))){
        std::cout<<"open file error\n";
    }

    string temp; //记录每一行的数据
    std::vector<std::string> v; 

    while(temp!="$Nodes"){
        getline(infile, temp, '\n');
    }
    getline(infile, temp, '\n');

    //输入node的总数
    outfile_node<<"\t"<<geo_v[3]<<"\n";
    //输入每一条node
    while(temp!="$EndNodes"){
        Split_String(temp,v);
        for(std::vector<string>::size_type i = 0; i != v.size(); ++i)
            outfile_node<<"\t"<<v[i];
        outfile_node<<"\n";
        getline(infile, temp, '\n');
        v.clear();
    }

    outfile_node.close();

}
/*int main()
{
    // 创建文件夹
    string folderpath = "D:\\Bridge Project\\ini\\Project\\swust";
    if (_access(folderpath.c_str(), 0)!=0)
        _mkdir(folderpath.c_str());
    
    //生成文件
    //复写不需改动的文件.txt
    string C_FGP = "D:\\Bridge Project\\swust\\CRWG_CSWG_VSIE_FGP.txt";
    string C_FGP_copy = "D:\\Bridge Project\\ini\\Project\\swust\\CRWG_CSWG_VSIE_FGP.txt";
    copyfile(C_FGP.c_str(),C_FGP_copy.c_str());

    string V_FGP = "D:\\Bridge Project\\swust\\CRWG_CSWG_VSIE_VGP.txt";
    string V_FGP_copy = "D:\\Bridge Project\\ini\\Project\\swust\\CRWG_CSWG_VSIE_VGP.txt";
    copyfile(V_FGP.c_str(),V_FGP_copy.c_str());

    //复写需要重写的文件.dat
    //读取ini文件信息
    string ini_file = "D:\\Bridge Project\\Detool\\mom.ini";
    INIReader reader(ini_file.c_str());
    
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<<ini_file<<"\n";
        return 1;
    }

    // //写入文件CRWG_CSWG_VSIE_Tet.DAT
    string Ctet_path = "D:\\Bridge Project\\ini\\Project\\swust\\CRWG_CSWG_VSIE_Tet.DAT";
    Ctet_file(reader,Ctet_path.c_str());

    //写入文件INCF_BiRCS.DAT
    string INCF_path = "D:\\Bridge Project\\ini\\Project\\swust\\INCF_BiRCS.DAT";
    Incf_file(reader, INCF_path.c_str());

    //写入文件WFB_BiRCS.DAT
    string Wfb_path = "D:\\Bridge Project\\ini\\Project\\swust\\WFB_BiRCS.DAT";
    Wfb_file(reader,Wfb_path.c_str());


    return 0;
}*/

// int main(){
//     // //打开要读取的文件
//     // string filename="D:\\Bridge Project\\1De_sw\\test.dat";
//     // string outfile="D:\\Bridge Project\\ini\\Project\\swust";//生成文件的目录
//     // std::ifstream infile ;

//     // infile.open(filename.c_str(),std::ios::out);
//     // if(!infile){
//     //     std::cout<<"open file error\n";
//     //     return 1;
//     // }

//     // string temp;
//     // getline(infile, temp, '\n');
//     // while(temp!="$Geo"){
//     //     getline(infile, temp, '\n');
//     // }

//     // //提取elem和node的数量存在geo_v中
//     // std::vector<std::string> geo_v;
//     // getline(infile, temp, '\n');
//     // while(temp!="$EndGeo"){
//     //     Split_String(temp,geo_v);
//     //     getline(infile, temp, '\n');
//     // }
    
//     // for(std::vector<string>::size_type i = 0; i != geo_v.size(); ++i)
//     //         std::cout<<geo_v[i]<<"\t";
//     // //生成node和elem文件
//     // //node_geo(infile, geo_v, outfile);
//     // elem_geo(infile, outfile);

//     // infile.close();

//     std::map<string,int> pre_;
//     // std::vector<string> key_pre = {"diagonal","normalized BF","no"};
//     // std::cout<<"size  "<<key_pre.size()<<"\n";

//     // for(int i=0; i<key_pre.size();i++){
//     //     pre_.insert(std::pair<string,int>(key_pre[i],i+1));
//     // }

//     // string test = "diagonal";
//     // if(pre_.find(test)==pre_.end())
//     //     std::cout<<"not find\n";

//     return 0;
// }