#include <io.h>
#include <direct.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>

#include "INIReader.h"
#include "INI_Handler.h"
#include <fstream>

#define M_PI 3.14159265358979323846

using std::string;

//字符串以什么开头和结尾
namespace std {
template <typename charT>
inline bool starts_with(const basic_string<charT>& big, const basic_string<charT>& small) {
    if (&big == &small) return true;
    const typename basic_string<charT>::size_type big_size = big.size();
    const typename basic_string<charT>::size_type small_size = small.size();
    const bool valid_ = (big_size >= small_size);
    const bool starts_with_ = (big.compare(0, small_size, small) == 0);
    return valid_ and starts_with_;
}

template <typename charT>
inline bool ends_with(const basic_string<charT>& big, const basic_string<charT>& small) {
    if (&big == &small) return true;
    const typename basic_string<charT>::size_type big_size = big.size();
    const typename basic_string<charT>::size_type small_size = small.size();
    const bool valid_ = (big_size >= small_size);
    const bool ends_with_ = (big.compare(big_size - small_size, small_size, small) == 0);
    return valid_ and ends_with_;
}
}  // namespace std

INI_Handler::INI_Handler(const string& filename_ini, const string& filename_out, const string& name_pos){
    ini_file = filename_ini;
    out_file = filename_out;
    name = name_pos;
    _error = 0;
     // 读mat文件
    _error = mat_read();

}

INI_Handler::~INI_Handler()
{
}

//复制不需要改动的文件
int INI_Handler::copyfile(const string& read_flie,string write_file){
    string temp;
    std::ifstream infile;
    std::ofstream outfile;

    infile.open(read_flie,std::ios::in);
    if(!infile)
        return -1;
    outfile.open(write_file,std::ios::out);
    if(!outfile)
        return -1;
    while (!infile.eof()){
        getline(infile, temp, '\n');
        outfile <<temp<<"\n";
    }
    outfile.close();
    infile.close();

    return 0;
}
int INI_Handler::C_FGP(){
    string file_fgp = ini_file + "CRWG_CSWG_VSIE_FGP.txt";
    std::ofstream odat_CTet;
    odat_CTet.open(file_fgp.c_str(),std::ios::out);

    odat_CTet<<"7   #for Volt and Ar_gauss_(NumEdg)=NSGau  1/3/4/6/7/9/12/16/25/36"<<"\n";
    odat_CTet<<"7     9    #PEC ,PSPS_no_singlarity=>NSGauO,PSPS_outer_singlarity=>NSGaui"<<"\n";
    odat_CTet<<"5     5        # Linsso,Linssi   !PEC面积分的奇异性处理情况下内外层线积分点"<<"\n";
    odat_CTet<<"5      #  Np_1D_ST    !面面重合线高斯点  #  处理PEC三角形单元近奇异性时线高斯选择*"<<"\n";
    odat_CTet<<"5    #  Np_1D_EA      !相邻三角形线高斯点"<<"\n";
    odat_CTet<<"5      #  Np_1D_VA    !共点三角形线高斯点"<<"\n";

    odat_CTet.close();
    return 0;
}
int INI_Handler::C_VGP(){
    string file_fgp = ini_file + "CRWG_CSWG_VSIE_VGP.txt";
    std::ofstream odat_CTet;
    odat_CTet.open(file_fgp.c_str(),std::ios::out);

    odat_CTet<<"5  8   #NdGaVvo,NdGaVvi  !描述电介质的体积分中体体积分外内层体高斯积分节点(选择4、5、8、10、11、15、27、35）"<<"\n";
    odat_CTet<<"5   5  5  # Linvv1,Linvv2,Linvv3 !体积分体体积分奇异性处理内3层线高斯积分节点"<<"\n";
    odat_CTet<<"5   5  5   #  Linvs1,Linvs2,Linvs3   !体面积分奇异性处理线高斯积分节点"<<"\n";
    odat_CTet<<"5   5  5    # NGLino1,NGLino2,NGLino3 !计算电通量密度矢量的线高斯积分节点"<<"\n";
    odat_CTet<<"7   7     # NGaVss_So  NGaVss_Si    !计算体积分中面面积分的外层面高斯积分节点(选择 1/3/4/6/7/9/12/16/25/36"<<"\n";
    odat_CTet<<"5   5    # Linvv_SS1,Linvv_SS2            !体积分中的边界处面面积分内层线高斯点"<<"\n";

    odat_CTet.close();
    return 0;
}
//找到当前文件夹下以目标后缀pos结尾的一个文件
void find_file_pos(string& filename, const string pos){
    //获取当前路径下的第一个mat文件
    string mat_file; 
	char   buffer[80];   
	getcwd(buffer, 80);
    string path = buffer;
    string p;
     //文件句柄  
    long   hFile   =   0;  
    //文件信息
    struct _finddata_t fileinfo;//文件各种信息的结构体  
    if((hFile = _findfirst(p.assign(path).append("\\*").append(pos).c_str(),&fileinfo)) !=  -1)  
    {  
        do  
        {  
            mat_file = p.assign("").append(fileinfo.name);
            break;

        }while(_findnext(hFile, &fileinfo)  == 0);  
        _findclose(hFile);  
    } 
    filename = mat_file;
}
//以空格分割字符串，不管连续空格的数量
void INI_Handler::Split_String(const std::string& s, std::vector<std::string>& v){
    string temp="";
    for(int i=0;i<s.length();i++){
        if((s[i]==' ')&&(temp.length()>0)){
            v.push_back(temp);
            temp = "";
        }
        else if(s[i]!=' '){
            temp += s[i];
        }
    }
    if(temp.size())
        v.push_back(temp);
}
void Split_String_elem(const std::string& s, std::vector<int>& v){
    string temp="";
    for(int i=0;i<s.length();i++){
        if((s[i]==' ')&&(temp.length()>0)){
            v.push_back(atof(temp.c_str()));
            temp = "";
        }
        else if(s[i]!=' '){
            temp += s[i];
        }
    }
    if(temp.size())
        v.push_back(atof(temp.c_str()));
    for(int i=v.size();i<12;i++)
        v.push_back(0);
    v.push_back(atof(temp.c_str()));
}
void Split_String_boundaries(const std::string& s, std::vector<int>& v){
    string temp="";
    for(int i=0;i<s.length();i++){
        if((s[i]==' ')&&(temp.length()>0)){
            v.push_back(atof(temp.c_str()));
            temp = "";
        }
        else if(s[i]!=' '){
            temp += s[i];
        }
    }
    if(temp.size())
        v.push_back(atof(temp.c_str()));
}
//以特定字符，切割字符串
void INI_Handler::SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c){
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
//生成CRWG_CSWG_VSIE_Tet.DAT文件/ini
int INI_Handler::Ctet_file(INIReader reader){
    string Ctet_path = out_file + "CRWG_CSWG_VSIE_Tet.DAT";
    // //判断文件是否存在
    // if(_access(Ctet_path.c_str(),0)==0)
    //     return 0;
       
    std::ofstream odat_CTet;
    odat_CTet.open(Ctet_path.c_str(),std::ios::out);
    if(!odat_CTet)
        return -2;
    odat_CTet<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<reader.GetReal("BACKGROUND_MEDIA", "EpsilonRe", 0)<<"    "
             <<reader.GetReal("BACKGROUND_MEDIA", "EpsilonIm", 0)<<"\t\t"<<'#'<<" BACKGROUND_MEDIA: EpsilonRe and EpsilonIm"<<"\n";
    odat_CTet<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<reader.GetReal("BACKGROUND_MEDIA", "MuRe", 0)<<"    "
             <<reader.GetReal("BACKGROUND_MEDIA", "MuIm", 0)<<"\t\t"<<'#'<<" BACKGROUND_MEDIA: MuRe and MuIm"<<"\n";

    odat_CTet.close();
    return 0;
}
int INI_Handler::Ctet_file_matr(std::map<string,int> matr_num,int tetrahedron_num){
    string Ctet_path = out_file + "CRWG_CSWG_VSIE_Tet.DAT";
    std::ofstream odat_CTet;
    odat_CTet.open(Ctet_path.c_str(),std::ios::app);
    if(!odat_CTet)
        return -2;
    odat_CTet<<"\n";

    //磁材料如果四面体数量为0，则加入0\n 0 0
    if(tetrahedron_num!=0){
       std::map<string,int>::iterator iter;
       int bhao_start = 1;
       int bhao_end = 0;
       int xuhao = 0;
       for(iter=matr_num.begin();iter!=matr_num.end();iter++){
           if(matr_num[iter->first]!=0){
                bhao_end += matr_num[iter->first];
                odat_CTet<<"\n"<<++xuhao<<"\n";
                odat_CTet<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<bhao_start<<"   "<<bhao_end<<"\n";
                odat_CTet<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<mat_value[iter->first]<<"   0.0"<<"\n";
                bhao_start += matr_num[iter->first];
           }
       } 
    }
    else{
        odat_CTet<<"0"<<"\n"<<"0   0"<<"\n";
    }
    //电材料
    odat_CTet<<"\n0"<<"\n"<<"0   0"<<"\n";

    odat_CTet.close();
    return 0;
}
//生成INCF_BiRCS.DAT文件/ini
//计算INCF传播方向的球坐标
void setPolizationDirection(double * inc_dir, double angle_, double * pol_dir){
        double r_dir[3];
        double z_dir[3] = {0, 0, 1.0};
        double theta_xyz[3];
        double phi_xyz[3];

        double theta = inc_dir[0] / 180.0 * M_PI;
        double phi = inc_dir[1] / 180.0 * M_PI;
        double angle =  angle_ / 180.0 * M_PI;
        
        r_dir[0] = std::sin(theta) * std::cos(phi);
        r_dir[1] = std::sin(theta) * std::sin(phi);
        r_dir[2] = std::cos(theta);

        CROSS_PRODUCT(z_dir, r_dir, phi_xyz);
        CROSS_PRODUCT(phi_xyz, r_dir, theta_xyz);

        for(int i = 0; i < 3; i++){        
            pol_dir[i] = theta_xyz[i] * std::cos(angle) 
                         - phi_xyz[i] * std::sin(angle);
        }

       // CROSS_PRODUCT(r_dir, e_pol_dir, h_pol_dir);
    }
int INI_Handler::Incf_file(INIReader reader,const string& incf_flie){
    double inc_dir[2],pol_dir[3];
    std::vector<std::string> v_theta;
    std::vector<std::string> v_phi;
    double angle;
    string split_=":";
    std::ofstream odat_INCF;

    odat_INCF.open(incf_flie,std::ios::out);
    if(!odat_INCF)
        return -3;
    
    //写入频率
    SplitString(reader.GetString("INCIDENT_WAVE", "Theta", "0:0:0"),v_theta,split_);
    SplitString(reader.GetString("INCIDENT_WAVE", "Phi", "0:0:0"),v_phi,split_);
    angle = reader.GetInteger("INCIDENT_WAVE", "PolarizationAngle", 0)*1.0;

    inc_dir[0] = atof(v_theta[0].c_str());
    inc_dir[1] = atof(v_phi[0].c_str());
    odat_INCF<<"    "<<"1.000"<<"\t\t\t"<<'#'
             <<"  Electric field amplitude of incident plane electromagnetic wave"<<"\n";
    //分别写入极化方向球坐标
    odat_INCF<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<"    "<<inc_dir[0]<<"    "
             <<inc_dir[1]<<"\t\t"<<'#'<<"  Spherical coordinates of polarization direction /° TheE0,PhiE0"<<"\n";
    //计算传播方向球坐标，并写入
    setPolizationDirection(inc_dir,angle,pol_dir);
    odat_INCF<<std::setiosflags(std::ios::fixed) <<std::setprecision(1)<<"    "<<pol_dir[0]<<"    "
             <<pol_dir[1]<<"\t\t"<<'#'<<" Spherical coordinates of propagation direction /°\n";
    
    odat_INCF.close();

    return 0;
}
//生成WFB_BiRCS.DAT文件/ini  
int INI_Handler::Wfb_file(INIReader reader,const string& wfb_flie,std::map<string,int>& pre_, std::map<string,int>& type_){
    std::vector<std::string> v_wfb;
    string split_=":";
    std::ofstream odat_Wfb;
    odat_Wfb.open(wfb_flie,std::ios::out);
    if(!odat_Wfb)
        return -4;
    //FREQUENCY
    int freq_num=0;
    SplitString(reader.GetString("FREQUENCY", "Value", "0:0:0"),v_wfb,split_);
    if(atof(v_wfb[1].c_str())!=0)
        freq_num = (atof(v_wfb[2].c_str())-atof(v_wfb[0].c_str()))/atof(v_wfb[1].c_str());

    odat_Wfb<<0.5<<"\t\t"<<'#'<<"  Integral weight of CFIE [1.0]=EFIE;  [0.0]=MFIE\n";
    odat_Wfb<<atof(v_wfb[0].c_str())<<"\t\t"<<'#'<<" Freq0 Initial frequency of incident wave\n";
    odat_Wfb<<atof(v_wfb[1].c_str())<<"\t\t"<<'#'<<" Dfreq Calculate the interval between frequencies\n";
    odat_Wfb<<freq_num<<"\t\t"<<'#'<<" Numfreq Calculate the numeber of frequencies\n";

    //SLOVER-Preconditioner
    string f_key; //关键字
    f_key = reader.GetString("SOLVER", "Preconditioner", "no");
    if(pre_.find(f_key)==pre_.end())
        f_key="no";
    odat_Wfb<<pre_[f_key]<<"    	# preconditioner：1->diagonal, 2->normalized BF, 3->no     #normalized BF produces very good results\n";
    
    //SLOVER-type
    f_key = reader.GetString("SOLVER", "Type", "LU");
    if(type_.find(f_key)==type_.end())
        f_key="LU";
    odat_Wfb<<type_[f_key]<<"    	# 1->Gauss, 2->CGM, 3->BCGM, 4->GMRES, 5->LU;\
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
    else {
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

    return 0;
}
//生成node.geo文件/geo
int INI_Handler::node_geo(std::ifstream& infile, std::vector<std::string>& geo_v, const string& out_file){
    string out_filename = out_file+ "node.geo";
    std::ofstream outfile_node ;

    outfile_node.open(out_filename.c_str(),std::ios::out);
    if((!outfile_node)||(!infile)||(geo_v.size()<4))
        return -5;

    string temp; //记录每一行的数据
    while(!infile.eof()){
        getline(infile, temp, '\n');
        if(temp=="$Nodes")
            break;
    }
    //输入node的总数
    outfile_node<<"  "<<geo_v[3]<<"\n";
    int num_of_nodes = 0;
    num_of_nodes = atoi(geo_v[3].c_str());
    //输入每一条node
    getline(infile, temp, '\n');
    for(int i=0;i<num_of_nodes;i++){
        outfile_node<<"  "<<temp<<"\n";
        getline(infile, temp, '\n');
    }
    outfile_node.close();

    return 0;
}
//生成elem.geo文件/geo
//读取mat文件
int INI_Handler::mat_read(){
    //获取当前路径下的第一个mat文件,//打开mat文件
    //string mat_pos = ".mat";
    string mat_name;
    mat_name = ini_file + name + ".mat";
    //find_file_pos(mat_name, mat_pos);
    std::ifstream mat_in;
    mat_in.open(mat_name.c_str(),std::ios::in);

    if(!mat_in)
        return -7;

    string temp;
    string mat_start, mat_end;
    string mat_num;
    double _value;
    std::vector<std::string> mat_v; 
    
    mat_start = "NUMBER:";
    mat_end = "VALUE";

    while(!mat_in.eof()){
        getline(mat_in,temp,'\n');
        if(std::starts_with(temp,mat_start)){
            Split_String(temp, mat_v);
            mat_num = mat_v[1];
        }
        if(std::starts_with(temp, mat_end)){
            Split_String(temp,mat_v);
            _value = atof(mat_v[1].c_str());
            if(!mat_num.empty())
               mat_value[mat_num]=_value; 
        }
        mat_v.clear();
    }
    return 0;
}
bool cmp_elem(const std::vector<int>& a, const std::vector<int>& b)
{
    return a[12] < b[12];
}
int INI_Handler::mid_elem(std::ifstream& infile, std::vector<std::string>& geo_v, const string& out_file){
    //count the number of material
    std::ofstream outfile_midelem ;
    outfile_midelem.open(out_file,std::ios::out);

     if((!outfile_midelem)||(!infile))
        return -6;
    //读取mat_value的数据
    std::map<string,double>::iterator iter;
    //材料的种类
    std::map<string,int> matr_num;
    //读取材料属性的种类，并初始化数量为0
    for(iter=mat_value.begin();iter!=mat_value.end();iter++)
        matr_num[iter->first] = 0;
    string temp; //记录每一行的数据

    while(!infile.eof()){
        getline(infile, temp, '\n');
        if(temp=="$Elements")
         break;
    }
    getline(infile, temp, '\n');
    getline(infile, temp, '\n');

     // find  Element
    std::vector<std::vector<int>> elem_4;
    std::vector<std::vector<int>> elem_2;
    std::vector<int> v; //读取需要写入的数据，做切割存到v中
    int number_of_elements = atoi(geo_v[1].c_str()); //elements总数
    for(int i = 0; i < number_of_elements; i++){
        Split_String_elem(temp,v);
        if(v[1]==4)
            elem_4.push_back(v);
        else
            elem_2.push_back(v);
        v.clear();
        getline(infile, temp, '\n');
    }

    //find $Boundaries
    while(!infile.eof()){
        getline(infile, temp, '\n');
        if(temp=="$Boundaries")
            break;
    }
    //找到boundaries数据的第一行，将boundaries加入到elem_2中
    getline(infile, temp, '\n');
    getline(infile, temp, '\n');

    int num_of_boudaries = 0;
    num_of_boudaries = atoi(temp.c_str());

    //boundaries的起始编号
    int numbers_of_elements_2;
    numbers_of_elements_2 = elem_2.size();

    std::vector<std::vector<int>> boudaries;
    v.clear();
    boudaries.resize(num_of_boudaries);
    for(int i = 0; i < num_of_boudaries; i++)
    {
        getline(infile, temp, '\n');
        Split_String_boundaries(temp,v);
        int boundaries_size=0;
        boundaries_size=v.size();
        if((boundaries_size>3)&&(boundaries_size < 10)){
            boudaries[i].push_back(2);
            boudaries[i].push_back(numbers_of_elements_2++);
            for(int j=boundaries_size-3;j<boundaries_size;j++)
                boudaries[i].push_back(v[j]);
            for(int q=5;q<12;q++)
                boudaries[i].push_back(0);
        }
        else{
            boudaries[i].push_back(2);
            boudaries[i].push_back(numbers_of_elements_2++);
            for(int j=boundaries_size-6;j<boundaries_size;j++){
                boudaries[i].push_back(v[j]);
            }
            for(int q=8;q<12;q++)
                boudaries[i].push_back(0);
        }
        v.clear();
    } 
    //写elem.geo文件
    //将out_elem.geo的内容复写到elem.geo中，并加入四面体和三角形的数量
    outfile_midelem <<"  "<<elem_4.size()<<"\n";
    outfile_midelem <<"  "<<elem_2.size()+boudaries.size()<<"\n";
    string material_bh;
    //排序
    stable_sort(elem_4.begin(), elem_4.end(), cmp_elem);
    //写四面体各属性的elements
    for(int i=0;i<elem_4.size();i++){
        material_bh = std::to_string(elem_4[i][12]) ;
        matr_num[material_bh]++;
        outfile_midelem<<"  "<<material_bh;
        outfile_midelem<<"  "<<i+1;
        for(int j=2;j<11;j++)
            outfile_midelem<<"  "<<elem_4[i][j];
        outfile_midelem<<"\n"<<"  "<<elem_4[i][11]<<"\n";
    }
    //写三角形各属性的elements
    for(int i=0;i<elem_2.size();i++){
        outfile_midelem<<"  "<<2;
        outfile_midelem<<"  "<<i+1;
        for(int j=2;j<11;j++)
            outfile_midelem<<"  "<<elem_2[i][j];
        outfile_midelem<<"\n"<<"  "<<elem_2[i][11]<<"\n";
    }
     //写boundaries各属性的elements
    numbers_of_elements_2 = elem_2.size();
    for(int i=0;i<boudaries.size();i++){
        outfile_midelem<<"  "<<2;
        outfile_midelem<<"  "<<++numbers_of_elements_2;
        for(int j=2;j<11;j++)
            outfile_midelem<<"  "<<boudaries[i][j];
        outfile_midelem<<"\n"<<"  "<<boudaries[i][11]<<"\n";
    }

    
    outfile_midelem.close();
    //补充CRWG_CSWG_VSIE_Tet.DAT的内容
    Ctet_file_matr(matr_num,elem_4.size());
    return 0;
}
int INI_Handler::elem_geo(std::ifstream& infile, std::vector<std::string>& geo_v, const string& out_file){
    string out_elem = out_file + "elem.geo";//最终生成的文件
    mid_elem(infile, geo_v, out_elem);
    return 0;
}
//生成所有ini文件
int INI_Handler::ini_out(){
    // 创建文件夹
    if (_access(out_file.c_str(), 0)!=0)
        _mkdir(out_file.c_str());
    //生成文件,复写不需改动的文件.txt
    C_FGP();
    C_VGP();

    //复写需要重写的文件.dat
    //读取ini文件信息
    string inifile = ini_file+ name +".ini";
    INIReader reader(inifile.c_str());
    
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<<ini_file<<"\n";
        return 1;
    }

    // //写入文件CRWG_CSWG_VSIE_Tet.DAT,放入geo_out中
    _error = Ctet_file(reader);
    if(_error!=0)
        return _error;

    //写入文件INCF_BiRCS.DAT
    string INCF_path = out_file + "INCF_BiRCS.DAT";
    _error = Incf_file(reader, INCF_path.c_str());
    if(_error!=0)
        return _error;

    //写入文件WFB_BiRCS.DAT
    //slover中的preconditioner和type初始化
    //Preconditioner对应的关键字，Preconditioner对应的关键字赋值
    std::map<string,int> pre_;
    std::vector<string> key_pre = {"diagonal","normalized BF","no"};

    for(int i=0; i<key_pre.size();i++){
        pre_.insert(std::pair<string,int>(key_pre[i],i+1));
    }
    //type对应的关键字，type对应的关键字赋值,将其存到map容器中
    std::map<string,int> type_;
    std::vector<string> key_type = {"Gauss","CGM","BCGM","GMRES","LU"};

    for(int i=0; i<key_type.size();i++)
        type_.insert(std::pair<string,int>(key_type[i],i+1));
    //type=5没用，先改为4
    type_["LU"] = 4;
        
    string Wfb_path = out_file + "WFB_BiRCS.DAT";
    _error = Wfb_file(reader,Wfb_path.c_str(),pre_,type_);
    if(_error!=0)
        return _error;

    return 0;
}
//生成所有geo文件
int INI_Handler::geo_out(){
    //获取当前路径下的第一个dat文件,//打开dat文件
    // string dat_pos = ".dat";
    string filename;
    //find_file_pos(filename, dat_pos);
    filename = ini_file + name + ".dat";
   
    std::ifstream infile ;
    infile.open(filename.c_str(),std::ios::out);
    if(!infile)
        return 2;

    string temp;
    getline(infile, temp, '\n');
    while(!infile.eof()){
        getline(infile, temp, '\n');
        if(temp=="$Geo")
            break;
    }

    //提取elem和node的数量存在geo_v中
    std::vector<std::string> geo_v;
    getline(infile, temp, '\n');
    while(!infile.eof()){
        Split_String(temp,geo_v);
        getline(infile, temp, '\n');
        if(temp=="$EndGeo")
            break;
    }
    if(infile.eof())
        return -5;

    //生成node.geo文件
    _error = node_geo(infile, geo_v, out_file);
    if(_error!=0){
         return _error;
    }
    //生成elem.geo文件
     _error=elem_geo(infile, geo_v, out_file);
     if(_error!=0){
         return _error;
    }

    infile.close();

    return 0;
}
