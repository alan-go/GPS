//
// Created by root on 7/13/18.
//

#include "SerialData.h"
#include "GNSS.h"

int SerialData::ParaseGGA( char* gga){
//    printf("GGA: %s\n", gga);
    if(sendGGA){
        gnss->rtkManager.SentGGA(gga,strlen(gga));
    }
    auto split = []( char* str,char c,vector<string> &result){
        char temp[16],*p = str;
        int k=0,n=0;
        while ((*p++)&&14!=n){
            if(*p==c){
                k=0;n++;
                result.push_back(temp);
                memset(temp,0,16);
            } else{
                temp[k++]=*p;
            }
        }
    };
    vector<string> val;
    split(gga,',',val);

    GnssTime time;
    double tod=0.0,lat=0.0,lon=0.0,hdop=0.0,alt=0.0,msl=0.0,ep[6],tt;
    char ns='N',ew='E',ua=' ',um=' ';
    Vector3d pos;
    int i,solq=0,nrcv=0;


    for (i=1;i<val.size();i++) {
        switch (i-1) {
            case  0: tod =atof(val[i].c_str()); break; /* time in utc (hhmmss) */
            case  1: lat =atof(val[i].c_str()); break; /* latitude (ddmm.mmm) */
            case  2: ns  =*val[i].c_str();      break; /* N=north,S=south */
            case  3: lon =atof(val[i].c_str()); break; /* longitude (dddmm.mmm) */
            case  4: ew  =*val[i].c_str();      break; /* E=east,W=west */
            case  5: solq=atoi(val[i].c_str()); break; /* fix quality */
            case  6: nrcv=atoi(val[i].c_str()); break; /* # of satellite tracked */
            case  7: hdop=atof(val[i].c_str()); break; /* hdop */
            case  8: alt =atof(val[i].c_str()); break; /* altitude in msl */
            case  9: ua  =*val[i].c_str();      break; /* unit (M) */
            case 10: msl =atof(val[i].c_str()); break; /* height of geoid */
            case 11: um  =*val[i].c_str();      break; /* unit (M) */
        }
    }
    if ((ns!='N'&&ns!='S')||(ew!='E'&&ew!='W')) {
        printf("invalid nmea gpgga format,%c,%c\n",ns,ew);
        return 0;
    }
    int week = gnss->utcTime->tm_wday;
    int hour = floor(tod/10000.0);
    tod-=hour*10000.0;
    int min = floor(tod/100.0);
    tod-=min*100.0;
    double secoday = 3600.0*hour+60.0*min+tod;

//    if (sol->time.time==0.0) {
//        trace(2,"no date info for nmea gpgga\n");
//        return 0;
//    }
    pos(0)=(ns=='N'?1.0:-1.0)*dmm2deg(lat)*D2R;
    pos(1)=(ew=='E'?1.0:-1.0)*dmm2deg(lon)*D2R;
//    pos(2)=alt+msl;
    pos(2)=alt;
    //<6> GPS状态， 0初始化， 1单点定位， 2码差分， 3无效PPS， 4固定解， 5浮点解， 6正在估算 7，人工输入固定值， 8模拟模式， 9WAAS差分
    Solution sol(GnssTime(secoday),pos);
    solRaw.push_front(sol);
    printf("tod,tod %f,%f\n",tod,GnssTime(secoday).tod);
    sol.Show("###GGA###");
//    time2epoch(sol->time,ep);
//    septime(tod,ep+3,ep+4,ep+5);
//    time=utc2gpst(epoch2time(ep));
//    tt=timediff(time,sol->time);
//    if      (tt<-43200.0) sol->time=timeadd(time, 86400.0);
//    else if (tt> 43200.0) sol->time=timeadd(time,-86400.0);
//    else sol->time=time;
//    LLA2XYZ(pos,sol->rr);
//    sol->stat=0<=solq&&solq<=8?solq_nmea[solq]:SOLQ_NONE;
//    sol->ns=nrcv;
//
//    sol->type=0; /* postion type = xyz */
//
//    trace(5,"decode_nmeagga: %s rr=%.3f %.3f %.3f stat=%d ns=%d hdop=%.1f ua=%c um=%c\n",
//          time_str(sol->time,0),sol->rr[0],sol->rr[1],sol->rr[2],sol->stat,sol->ns,
//          hdop,ua,um);
    return 0;
}

SerialData::SerialData() :
flag(0),count(0),lengthNMEA(0),lengthUBX(0),baudRate(115200),stopCapture(false),showData(false){
    sp_ = nullptr;
}

SerialData::~SerialData() {

}

void SerialData::StartCapture(const std::string serialPort, unsigned int baudRate, char *saveName) {
    sprintf(saveName,"../data/device%d_%s.data",id,gnss->timeName);
    std::ofstream outF;
    if(logOpen)outF.open(saveName,std::ofstream::binary);
    try {
        boost::asio::io_service ios;
        sp_ = new boost::asio::serial_port(ios, serialPort);
        sp_->set_option ( boost::asio::serial_port::baud_rate ( baudRate ) );
        printf ( "successfully opened port %s\n", serialPort.c_str() );

        while ( !stopCapture ) {
            char tmp[256];
            auto transferred = sp_->read_some ( boost::asio::buffer ( tmp ) );
            if ( transferred <= 0 ) {
                usleep ( 100 );
                printf ( "serial port return 0\n" );
                continue;
            }
//            printf ( "transferred = %d , flag = %d\n",transferred, flag  );
            printf("get %d byte from USB%d  ::: %s\n",transferred, id,tmp);
            if(logOpen)outF.write(tmp,transferred);
            ScanSerialData(tmp,transferred);
        }
        sp_->close();
        outF.close();
    } catch ( ... ) {
        printf ( "failed to open serial port\n" );
        return;
    }
}

void SerialData::ScanSerialData(char *tmp, int transferred) {
    for ( int i = 0; i<transferred; i++ ) {
        if ( '$'==tmp[i]&&'G'==tmp[i+1] ) {
            flag = 1;
            lengthNMEA = 0;
            //printf ( "\n\nflag -> 1, i = %d\n",i );
        }
        if ( ( u_char ) tmp[i] == 0xB5 && ( u_char ) tmp[i+1] == 0x62 ) {
            flag = 2;
            lengthUBX = 0;
            //printf ( "\n\nflag -> 2, i = %d\n",i );
        }
        if(1==flag){
            bufferNMEA[lengthNMEA] = tmp[i];
            lengthNMEA++;
        }
        if(2==flag){
            bufferUBX[lengthUBX] = tmp[i];
            if(6==lengthUBX){
                //length is writen in 4,5 in protocol
                lengthUBXProtocol = *(u_int16_t*)(bufferUBX+4);
            }
            lengthUBX++;
        }
        if(('\n'==tmp[i])&&1==flag){
//            if(('\n'==tmp[i]||'\r'==tmp[i])&&1==flag){
            printf("Got a NMEA,l = %d.:%s",lengthNMEA,bufferNMEA);
            if(!memcmp(bufferNMEA+3,"GGA",3)&&lengthNMEA<128)ParaseGGA(bufferNMEA);
//            if (showData)
            if (1)
                for(int k = 0;k<lengthNMEA;k++){
                    printf("%c",bufferNMEA[k]);
                }
            flag=0;
        }
        //there are 8 extra bytes besides the playload;
        if(lengthUBX==lengthUBXProtocol+8 && 2==flag){
//            printf("\nGot a UBX %02x %02x,l = %d, i = %d, count = %d\n",bufferUBX[2],
//                   bufferUBX[3], lengthUBXProtocol, i,count++);
            if(showData && lengthUBX<1024)
                for(int k = 0;k<lengthUBX;k++){
                    printf("%02x ",(u_char) bufferUBX[k]);
                }
            parse_UBX(bufferUBX);
            flag=0;
        }
    }
}

void SerialData::parse_UBX(char *buffer) {
    if(0x02==(u_char)buffer[2]){
        if(1!=mode) return;
        if(0x15==(u_char)buffer[3]){
//            printf("\n--0--ParseRawData,len = %d\n",lengthUBX);
            gnss->ParseRawData(buffer, lengthUBX);
//            printf("\n--1--ParseRawData,%d\n",gnss->useBeiDou);
        }
        if(0x13==(u_char)buffer[3]){
//            printf("\n----0----ParseBstSubFrame,%d\n",gnss->useBeiDou);
            gnss->svsManager.UpdateEphemeris(buffer);
//            printf("\n----1----ParseBstSubFrame,%d\n",gnss->useBeiDou);
        }
    }

    if(0x01==(u_char)buffer[2]){
        if(0x06==(u_char)buffer[3]){
            if(lengthUBX>64)
                return;
           parse_UBX_0106(buffer);
        }
    }
}

void SerialData::parse_UBX_0106(char *buffer) {
    char* data = buffer+6;
    double tu[Nsys]={0};
    uint32_t itow = *(uint32_t*)data;
    int32_t ftow = *(int32_t*)(data+4);
    int16_t week = *(int16_t *)(data+8);
    GnssTime time(int(week),double(itow*1e-3+ftow*1e-9));
    Vector3d xyz(*(int32_t*)(data+12),*(int32_t*)(data+16),*(int32_t*)(data+20));
    Vector3d vxyz(*(int32_t*)(data+28),*(int32_t*)(data+32),*(int32_t*)(data+36));
    Solution sol(time,xyz,vxyz,tu);
    solRaw.push_front(sol);
}

int SerialData::StopDevice() {
    stopCapture = true;
    try {
        if (thread > 0) {//(void*)0
            pthread_join(thread, nullptr);
        }
        thread = 0;
        if (sp_ != nullptr) {
            sp_->close();
            delete sp_;
        }
        sp_ = nullptr;
    }catch (...){}

}

int SerialData::WtiteSerial() {
    return 0;
}
