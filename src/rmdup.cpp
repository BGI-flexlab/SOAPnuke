//
// Created by berry on 2020-10-20.
//


#include "rmdup.h"


rmdup::rmdup(uint64_t* data, uint64_t size) {
    this->data=data;
    this->size=size;
}

void rmdup::markDup(bool *dupFlag) {
    uint32_t prime=getPrime();
//    cout<<get_local_time()<<"\tprime:\t"<<prime<<endl;
    uint32_t *count=new uint32_t[prime];
    memset(count,0,sizeof(uint32_t)*prime);
//    unordered_set<string> tmpset;
    uint32_t contentSize=0;
    uint32_t multisize=0;
//    vector<uint64_t> dupIdx,dupIdx2;
//    unordered_set<uint64_t> tmpset2,dupset;
    for(uint64_t i=0;i<size;i++)
    {
//        string curStr=md52string(data[i]);
//        if(tmpset.find(curStr)!=tmpset.end()){
//            checkDupNum++;
////            dupset.insert(curStr);
//            cout<<"here\t"<<curStr<<endl;
//            dupIdx.push_back(i);
//        }else {
//            tmpset.insert(curStr);
//        }
        uint32_t idx=data[i]%prime;
        count[idx]++;
    }
//    cout<<get_local_time()<<"\tfirst scanning done"<<endl;
//    cout<<checkDupNum<<endl;
    uint32_t *multiIdx=new uint32_t[prime];
    memset(multiIdx,0,sizeof(uint32_t)*prime);
    uint32_t tmpIter=0;
    for(uint32_t i=0;i<prime;i++)
    {
        if(count[i]>1)
        {
            multiIdx[i]=tmpIter;
            tmpIter++;
            contentSize++;
            multisize+=count[i];
        }
    }
//    cout<<get_local_time()<<"\tsecond scanning done"<<endl;
    uint64_t **dataContent=new uint64_t *[contentSize];
//    memset(dataContent,0,sizeof(uint64_t*)*contentSize);
    uint32_t contentIdx=0;
    for(uint32_t i=0;i<prime;i++)
    {
        if(count[i]>1)
        {
            dataContent[contentIdx]=new uint64_t[count[i]];
            contentIdx++;
        }
    }
//    cout<<get_local_time()<<"\tthird scanning done"<<endl;

    int *curCount=new int[contentSize];
    memset(curCount,0,sizeof(int)*contentSize);

    for(uint64_t i=0;i<size;i++)
    {
        uint32_t idx=data[i]%prime;
        if(i%10000000==0)
        {
//            cout<<get_local_time()<<"\tin forth scanning:\t"<<i<<endl;
        }
        if(count[idx]>1)
        {
            if(multiIdx[idx]>=contentSize)
            {
                cerr<<"Error:code error"<<endl;
                exit(1);
            }
            dataContent[multiIdx[idx]][curCount[multiIdx[idx]]]=data[i];
            curCount[multiIdx[idx]]++;
            if(curCount[multiIdx[idx]]==count[idx])
            {
                unordered_set<uint64_t> exists;
//                string curStr;
                for(uint32_t j=0;j<count[idx];j++)
                {
                    uint64_t curDa=dataContent[multiIdx[idx]][j];
//                    curStr=md52string(dataContent[multiIdx[idx]][j]);
                    if(exists.find(curDa)==exists.end())
                    {
                        exists.insert(curDa);
                    }else
                    {
//                        out<<dataContent[multiIdx[idx]][j]<<endl;
                        dataContent[multiIdx[idx]][j]=-1;
                    }
                }
            }
        }
    }
//    cout<<get_local_time()<<"\tforth scanning done"<<endl;
    memset(curCount,0,sizeof(int)*contentSize);
    uint32_t iter=0;
    for(uint64_t i=0;i<size;i++) {
        uint32_t idx=data[i]%prime;
        if (count[idx] > 1) {
            if (multiIdx[idx] >= contentSize) {
                cerr << "Error:code error" << endl;
                exit(1);
            }
            if (dataContent[multiIdx[idx]][curCount[multiIdx[idx]]] == -1) {
                dupFlag[i]=true;
//                dupIdx2.push_back(i);
                iter++;
            }
            curCount[multiIdx[idx]]++;
        }
    }
//    cout<<get_local_time()<<"\tfifth scanning done"<<endl;
//    for(vector<uint64_t>::iterator ix=dupIdx.begin();ix!=dupIdx.end();ix++){
//        cout<<*ix<<endl;
//    }
    for(uint64_t i=0;i<size;i++)
    {
        uint32_t idx=data[i]%prime;
        if(count[idx]>1)
        {
            if(dataContent[multiIdx[idx]]!=NULL)
            {
//                for(uint32_t j=0;j<count[idx];j++){
//                    delete[] dataContent[multiIdx[idx]][j];
//                    dataContent[multiIdx[idx]][j]=NULL;
//                }
                delete[] dataContent[multiIdx[idx]];
                dataContent[multiIdx[idx]]=NULL;
            }
        }
    }
//    cout<<get_local_time()<<"\tsixth scanning done"<<endl;
    delete[] dataContent;
    delete[] multiIdx;
    delete[] curCount;
    delete[] count;
}

uint32_t rmdup::getPrime(){
    uint32_t realSize=size>pow(2,32)-1?pow(2,32)-1:size;
    uint32_t prime=0;
    if(size>0&&size<10)
    {
        return size;
    }
//    float low=(float)49/50;
//    float high=1;
    uint32_t curNum=realSize;
    while(curNum--)
    {
        bool isPrime=true;
        for(uint32_t j=2;j<=sqrt(curNum);j++)
        {
            if(curNum%j==0)
            {
                isPrime=false;
                break;
            }
        }
        if(isPrime)
        {
            prime=curNum;
            break;
        }
        if(prime>0)
        {
            break;
        }
    }
    if(prime<=0)
    {
        cerr<<"Error:code error, please contact with "<<AUTHOREMAIL<<endl;
        exit(1);
    }
    return prime;
}

rmdup::~rmdup() {
//    for(uint64_t i=0;i<size;i++){
//        delete[] data[i];
//        data[i]=NULL;
//    }
    delete[] data;
    data=NULL;
}

uint64_t rmdup::getlong(unsigned char *str){
    uint64_t value=0;
    for(uint32_t i=0;i<8;i++)
    {
        value+=(uint64_t)str[i]<<(8*(7-i));
    }
    return value;
}

string rmdup::md52string(unsigned char *str) {
    string value;
    unsigned char *tmpStr=str;
    value=to_string(getlong(str));
    tmpStr+=8;
    value+=to_string(getlong(tmpStr));
    return value;
}
