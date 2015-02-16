#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

struct CTransStruct{
    uint64_t trans_id;
    uint64_t value;
    CTransStruct(uint64_t tid, uint64_t v): trans_id(tid), value(v) {}
    inline bool operator< (const CTransStruct& rhs) const {
        /* do actual comparison */ 
        if (this->trans_id < rhs.trans_id) return true;
        else if (this->trans_id > rhs.trans_id) return false;
        else return this->value < rhs.value;
    }
    inline bool operator> (const CTransStruct& rhs) const {return rhs < *this;}
    inline bool operator<=(const CTransStruct& rhs) const {return !(*this > rhs);}
    inline bool operator>=(const CTransStruct& rhs) const {return !(*this < rhs);}
};
bool combCTRS_upper(uint64_t trans_id, const CTransStruct& a) {
    return trans_id <= a.trans_id;
}
bool combCTRS_lower(const CTransStruct& a, uint64_t trans_id) {
    return a.trans_id < trans_id;
}
struct CTRSLessThan_t {
    bool operator() (const CTransStruct& left, const CTransStruct& right) {
        return left.trans_id < right.trans_id;
        //if (left.trans_id < right.trans_id) return true;
        //else if (right.trans_id < left.trans_id) return false;
        //else return left.value < right.value;
    }
    bool operator() (const CTransStruct& o, uint64_t target) {
        return o.trans_id < target;
    }
    bool operator() (uint64_t target, const CTransStruct& o) {
        return target < o.trans_id;
    }
} CTRSLessThan;

int main() {
    cout << "testing" << endl;

    vector<CTransStruct> v;
    vector<pair<uint64_t, uint64_t> > vp;

    for(int i=0; i<10; i++) {
        v.push_back(CTransStruct(i, i*1000));
        v.push_back(CTransStruct(i, i*1000+1));
        vp.push_back(pair<uint64_t,uint64_t>(i, i*1000));
        vp.push_back(pair<uint64_t, uint64_t>(i, i*1000+1));
    }

    for(unsigned int i=0; i<v.size(); i++) {
        cout << "[" << v[i].trans_id << ":" << v[i].value << "] ";
    }
    cout << endl;
    for(unsigned int i=0; i<v.size(); i++) {
        cout << "[" << vp[i].first << ":" << vp[i].second << "] ";
    }
    cout << endl;

    CTransStruct o1(2, 0);
    CTransStruct o11(2, -1);
    CTransStruct o2(5, -1);

    auto lb = lower_bound(v.begin(), v.end(), o1, CTRSLessThan);
    auto ub = upper_bound(v.begin(), v.end(), o1, CTRSLessThan);
    auto ub1 = upper_bound(v.begin(), v.end(), o11, CTRSLessThan);
    auto ub2 = upper_bound(v.begin(), v.end(), o2, CTRSLessThan);

    cout << "lb: 2 " << (lb == v.end()) << " - " << lb->value << endl;
    cout << "ub1: 2 " << (ub1 == v.end()) << " - " <<  ub1->value << endl;
    cout << "ub: 2 " << (ub == v.end()) << " - " <<  ub->value << endl;
    cout << "ub: 5 " << (ub2 == v.end()) << " - " << ub2->value << endl;
    
    pair<uint64_t, uint64_t> p1(2, 0);
    pair<uint64_t, uint64_t> p2(5, -1);

    auto plb = lower_bound(vp.begin(), vp.end(), p1);
    auto pub = upper_bound(vp.begin(), vp.end(), p1);
    auto pub2 = upper_bound(vp.begin(), vp.end(), p2);

    cout << "lb: 2 " << (plb == vp.end()) << " - " << plb->second << endl;
    cout << "ub: 2 " << (pub == vp.end()) << " - " <<  pub->second << endl;
    cout << "ub: 5 " << (pub2 == vp.end()) << " - " << pub2->second << endl;

    return 0;
}


