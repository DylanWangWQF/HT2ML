#include <string.h>
#include <NTL/matrix.h>
#include <vector>
#include <chrono>

#include "dispatcher.h"
#include "trace.h"

//HTparams* params
ecall_dispatcher::ecall_dispatcher()
{
    // e_HTparams = params;
}

int ecall_dispatcher::enclave_init(uint8_t* hecontext, size_t context_len, uint8_t** rawCtxt, size_t* rawCtxt_len, size_t* num_rawCtxt)
{
    int ret = 0;
    // if (context_len == 0)
    // {
    //     TRACE_ENCLAVE("passed hecontext is null string, context_len =0");
    //     ret = 1;
    //     goto exit;
    // }
    auto start= chrono::steady_clock::now();

    TRACE_ENCLAVE("Enclave: receive the HE params");
    stringstream ess;
    string e_hecontext = string(hecontext, hecontext + context_len);
    // cout << "Enclave: check hecontext: " << hecontext << endl;
    // cout << "Enclave: check context_len, containing context, sk and pk: " << context_len << endl;
    ess << e_hecontext;
    TRACE_ENCLAVE("Enclave: readfrom the string and reconstruct the context: ");
    e_context = Context::readPtrFrom(ess);
    TRACE_ENCLAVE("Enclave: Setup Seckey by reading from ss!");
    activeSecKey = make_unique<SecKey>(SecKey::readFrom(ess, *e_context));
    TRACE_ENCLAVE("Enclave: Setup PubKey by reading from ss!");
    activePubKey = make_unique<PubKey>(PubKey::readFrom(ess, *e_context));

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Enclave: Total Setup Time for HE = " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;

    // Initialize HT params
    // allocate memory for struct
    TRACE_ENCLAVE("Enclave: allocate memory for struct.");
    // struct contains vector, cannot use malloc.
    e_HTparams = new HTparams();
    // e_HTparams = (struct _HTparams*)malloc(sizeof(struct _HTparams));
    // if (e_HTparams == nullptr)
    // {
    //     ret = 1;
    //     TRACE_ENCLAVE("allocate e_HTparams failed with %d", ret);
    //     goto exit;
    // }

    // Begin: hard-code dataset info
    TRACE_ENCLAVE("Enclave: Initialize HT params.");
    e_HTparams->numFeatures = 7;
    e_HTparams->numClasses = 2;
    e_HTparams->numBits = 16;
    e_HTparams->numAllCategories = {2, 2, 2, 2, 2, 2, 2, 2};
    // cout << "Enclave: Check numClasses:" << e_HTparams->numClasses << endl;

    unordered_map<size_t, string> umap;
    for (size_t i = 0; i < (e_HTparams->numFeatures + 1); i++)
    {
        umap.clear();
        umap.insert(pair<size_t, string>(0, "10"));
        umap.insert(pair<size_t, string>(1, "01"));
        e_HTparams->FeaCatMapping.push_back(umap);
    }

    e_HTparams->successProbability = 0.95;
    e_HTparams->maxSamples = 75;
    e_HTparams->checkInterval = 60;
    e_HTparams->minSamples = 60;
    // Begin: hard-code dataset info


    TRACE_ENCLAVE("Enclave: Initialize the leafnodes list!");
    HoeffdingTree* GiniHoeffdingTree = new HoeffdingTree(e_HTparams->numFeatures, 
                                                         e_HTparams->numClasses, 
                                                         e_HTparams->numBits, 
                                                         e_HTparams->numAllCategories, 
                                                         e_HTparams->FeaCatMapping, 
                                                         e_HTparams->successProbability, 
                                                         e_HTparams->maxSamples, 
                                                         e_HTparams->checkInterval, 
                                                         e_HTparams->minSamples);
    leafnodes.push_back(GiniHoeffdingTree);

    TRACE_ENCLAVE("Enclave: Generate Root Tree Matrix!");
    vector<string> Rootmatrix;
    GiniHoeffdingTree->GenerateCurrentTreeMatrix(Rootmatrix);

    TRACE_ENCLAVE("Enclave: Encrypt Root Tree Matrix!");
    vector<Ctxt> ReturnedRootCtxts;
    EncryptMat(ReturnedRootCtxts, Rootmatrix);

    TRACE_ENCLAVE("Enclave: copy memory to host and return the Root Ctxts to host!");
    // num_returnedMat: number of Mat returned to host
    num_returnedMat = ReturnedRootCtxts.size();

    ess.str(std::string());
    ess.clear();
    for (size_t i = 0; i < num_returnedMat; ++i)
    {
        ReturnedRootCtxts[i].writeTo(ess);
    }
    string rctxt = ess.str();
    uint8_t* host_buf = (uint8_t*) oe_host_malloc(size_t(rctxt.size() + 1));
    memcpy(host_buf, (uint8_t*)rctxt.c_str(), rctxt.size() + 1);
    *rawCtxt = host_buf;
    *rawCtxt_len = rctxt.size();
    *num_rawCtxt = num_returnedMat;

exit:
    TRACE_ENCLAVE("Enclave: free memory in enclave_init()!");
    e_hecontext.shrink_to_fit();
    rctxt.shrink_to_fit();
    ReturnedRootCtxts.clear();
    Rootmatrix.clear();
    vector<Ctxt>().swap(ReturnedRootCtxts);
    vector<string>().swap(Rootmatrix);
    return ret;
}

int ecall_dispatcher::multipleCtxtsTransform(uint8_t* ectxt, size_t ectxt_len, size_t num_ectxt, size_t batchIdx, uint8_t** octxt, size_t* octxt_len, size_t* num_octxt)
{
    stringstream css;
    vector<Ctxt> ReceivedCtxts(num_ectxt, Ctxt(*activePubKey));
    // ReceivedCtxts = vector<Ctxt>(num_hostCtxt, Ctxt(*activePubKey));
    // cannot confirm the size of ReturnedCtxts at the beginning
    vector<Ctxt> ReturnedCtxts;

    // uint8_t* hostCtxt_cpy[num_resCtxt];
    // copy_arr_to_enclave(hostCtxt_cpy, num_resCtxt, hostCtxt, hostCtxt_length);

    TRACE_ENCLAVE("Enclave: Receiving Ctxt start!");
    cout << "Enclave: Check num_ectxt: " << num_ectxt << endl;
    // cout << "Enclave: Check ectxt_len: " << ectxt_len << endl;
    // cout << "Enclave: Check ectxt pointer: " << ectxt << endl;
    string etemp = string(ectxt, ectxt + ectxt_len);
    css << etemp;
    for (size_t i = 0; i < num_ectxt; i++)
    {
        TRACE_ENCLAVE("Enclave: %lu-th readFrom Ctxt!", i);
        //Note, if push_back without set num at the beginning, it's null ctxt!!!!!!!!!
        // helib::Ctxt newCtxt = helib::Ctxt::readFrom(css, *activePubKey);
        // ReceivedCtxts[i] = newCtxt;
        ReceivedCtxts[i].Ctxt::read(css);

        // cout << "Decrypt ReceivedCtxts[" << i << "] in multipleCtxtsTransform(): " << endl;
        // vector<long> cmsg;
        // e_hecontext->ea.decrypt(ReceivedCtxts[i], *activeSecKey, cmsg);
        // long k = 0;
        // for(long i = 0; i < 256; ++i){
        //     cout << cmsg[i] << " ";
        // }
        // cout << endl;
    }

    // invoke the training procedure
    TRACE_ENCLAVE("Enclave: invoke HTtrain()!");
    HTtrain(ReceivedCtxts, ReturnedCtxts, batchIdx);

    TRACE_ENCLAVE("Enclave: copy memory to host and return the Ctxts to host!");
    num_returnedMat = ReturnedCtxts.size();
    cout << "Enclave: Check ReturnedCtxts.size(): " << num_returnedMat << endl;

    css.str(std::string());
    css.clear();
    for (size_t i = 0; i < num_returnedMat; ++i)
    {
        TRACE_ENCLAVE("Enclave: %lu-th writeTo Ctxt!", i);
        ReturnedCtxts[i].writeTo(css);
    }
    string otemp = css.str();
    // free host_buf??
    uint8_t* host_buf = (uint8_t*) oe_host_malloc(size_t(otemp.size() + 1));
    memcpy(host_buf, (uint8_t*)otemp.c_str(), otemp.size() + 1);
    *octxt = host_buf;
    *octxt_len = otemp.size();
    *num_octxt = num_returnedMat;

    return 0;
}

int ecall_dispatcher::HTtrain(vector<Ctxt>& ReceivedCtxts, vector<Ctxt>& ReturnedCtxts, size_t batchIdx)
{
    list<HoeffdingTree*>::iterator iter;
    // get the current tree matrix
    vector<string> Allleafmatrix;
    //countedMat->current leaf nodes in list
    //map<size_t, Mat<size_t>> -> in each leaf node, <feature, feature matrix>
    vector<map<size_t, Mat<size_t>>> countedMat;
    // leafnodes.size(), <countIndice, SizeMapping>
    // SizeMapping = number of cols for each leaf node in the resultmatrix, all of SizeMapping = Allleafmatrix.size()
    vector<pair<size_t, size_t>> indiceAndsize;
    // Setup the map to record: key->feature, value->Mat for each feature.
    map<size_t, Mat<size_t>> mapleaf;
    // record for the SizeMapping
    size_t stemp;
    Mat<size_t> mtemp;
    // locate the leadnode index in list
    size_t numleaf = 0;

    TRACE_ENCLAVE("Enclave: Initialize container start!");
    for (iter = leafnodes.begin(); iter != leafnodes.end(); iter++)
    {
        stemp = 0;
        if ((*iter)->remFeatures.empty())
        {
            TRACE_ENCLAVE("Enclave: Only genereate label matrix!");
            mtemp.SetDims(e_HTparams->numClasses, 1);
            MatInit(mtemp);
            // pair<labelIndex, Mat>, labelIndex = numFeatures
            mapleaf.insert(pair<size_t, Mat<size_t>>(e_HTparams->numFeatures, mtemp));
            mtemp.kill();
            stemp = e_HTparams->numClasses;
        }else{
            for(auto it = ((*iter)->remFeatures).begin(); it != ((*iter)->remFeatures).end(); it++){
                mtemp.SetDims(e_HTparams->numClasses, e_HTparams->numAllCategories[(*it)]);
                MatInit(mtemp);
                // <feature, feature matrix>, mat = sufficientStatistics
                mapleaf.insert(pair<size_t, Mat<size_t>>((*it), mtemp));
                mtemp.kill();
                // sum the number of feature string in current leafnode
                stemp += ((e_HTparams->numClasses) * e_HTparams->numAllCategories[(*it)]);
            }
        }
        // cout << "Current leafnode's mapleaf size (remain features) after generate matrix: " << mapleaf.size() << endl;
        // cout << "Current leafnode's stemp in resultMatrix Cols after generate matrix: " << stemp << endl;
        // one mapleaf -> one leafnode
        countedMat.push_back(mapleaf);
        mapleaf.clear();
        //CountIndice should be 2
        indiceAndsize.push_back(make_pair(((*iter)->CountIndice()), stemp));
    }
    TRACE_ENCLAVE("Enclave: Initialize container end!");

    TRACE_ENCLAVE("Enclave: Decrypt Ctxts and record statistics in container!");
    DecryptMat(ReceivedCtxts, countedMat, indiceAndsize, batchIdx);

    TRACE_ENCLAVE("Enclave: HT training start!");
    size_t leafnodesCount = leafnodes.size();
    cout << "const leafnodes size before train: " << leafnodes.size() << endl;
    // leafnodes size is changing due to erase
    for (iter = leafnodes.begin(); iter != leafnodes.end();)
    {
        if (leafnodesCount-- == 0)
        {
            break;
        }
        size_t isErase = 0;
        (*iter)->Train(countedMat[numleaf], leafnodes, isErase);
                
        if(isErase == 1)
        {
            // here, not iter++
            leafnodes.erase(iter++);
            cout << "Confirm leafnodes size after earse: " << leafnodes.size() << endl;
        }else{
            iter++;
        }
        numleaf++;
    }
    countedMat.clear();
    indiceAndsize.clear();
    TRACE_ENCLAVE("Enclave: HT training end!");

    for (iter = leafnodes.begin(); iter != leafnodes.end(); iter++)
    {
        TRACE_ENCLAVE("Enclave: Generate Matrix!");
        (*iter)->GenerateCurrentTreeMatrix(Allleafmatrix);
    }

    TRACE_ENCLAVE("Enclave: Encrypt the Mat!");
    EncryptMat(ReturnedCtxts, Allleafmatrix);

    return 0;
}

void ecall_dispatcher::DecryptMat(vector<Ctxt>& ReceivedCtxts, vector<map<size_t, Mat<size_t>>>& countedMat, vector<pair<size_t, size_t>>& indiceAndsize, size_t batchIdx)
{
    if (countedMat.size() != indiceAndsize.size())
    {
        cout << "Errors in Set up countedMat and indiceAndsize, both of should be leafnodes number" << endl;
        return;
    }
    Mat<size_t> resultMatrix;
    // num_HostCtxt * dim is set in host, num_HostCtxt is set in host (64 / nrows = 4/2/1)
    cout << "num_returnedMat in last return: " << num_returnedMat << endl;
    resultMatrix.SetDims(/*num_HostCtxt * dim*/(batchIdx * dim), num_returnedMat * dim);
    MatInit(resultMatrix);
    vector<long> cmsg;

    for (size_t n = 0; n < batchIdx; ++n)
    {
        size_t idx = 0;
        for (size_t k = (n * num_returnedMat); k < ((n + 1) * num_returnedMat); ++k)
        {
            cmsg.clear();
            e_context->getEA().decrypt(ReceivedCtxts[k], *activeSecKey, cmsg);

            size_t l = 0;
            for(size_t i = (n * dim); i < ((n + 1) * dim); ++i)
            {
                for(size_t j = (idx * dim); j < ((idx + 1) * dim); ++j)
                {
                    resultMatrix[i][j] = cmsg[l];
                    l++;
                }
            }
            idx++;
        }
    }

    // for (size_t k = 0; k < num_returnedMat; ++k)
    // {
    //     cmsg.clear();
    //     e_context->getEA().decrypt(ReceivedCtxts[k], *activeSecKey, cmsg);

    //     size_t l = 0;
    //     for(size_t i = 0; i < dim; ++i){
    //         for(size_t j = (k * dim); j < ((k + 1) * dim); ++j){
    //             // long to size_t, correct? (unsigned) cmsg[l]
    //             resultMatrix[i][j] = cmsg[l];
    //             l++;
    //         }
    //     }
    // }
    
    queue<size_t> desiredCount;
    size_t startIndex;
    size_t endIndex;
    size_t recordIndex = 0;
    for (size_t k = 0; k < indiceAndsize.size(); ++k){
        // cout << k << "-th leafnode's <count, number>: " << "<" << indiceAndsize[k].first << ", " << indiceAndsize[k].second << ">" << endl;
        //! Record the desired count. The size of Queue = countedMat.NumCols(); indiceAndsize[k].second=k-th leafnode's feature string
        //calculate index
        startIndex = recordIndex;
        endIndex = recordIndex + indiceAndsize[k].second;
        
        for (size_t i = startIndex; i < endIndex; ++i) {
            size_t ctemp = 0;
            // resultMatrix.NumRows()= number of data samples
            for (size_t j = 0; j < (unsigned)resultMatrix.NumRows(); ++j) {
                if (resultMatrix[j][i] == indiceAndsize[k].first) {
                    ctemp++;
                }
            }
            desiredCount.push(ctemp);
        }
        recordIndex = endIndex;
    }

    // Assign the first value in the queue to the feature matrix in each leaf node.
    map<size_t, Mat<size_t>> ::iterator mapiter;
    for (size_t k = 0; k < countedMat.size(); k++){
    // k-th leafnode
        for (mapiter = countedMat[k].begin(); mapiter != countedMat[k].end(); mapiter++){
            // NumCols() = categories (sunny, rainy, windy); NumRows() = classes
            for (size_t i = 0; i < (unsigned)mapiter->second.NumCols(); i++) {
                for (size_t j = 0; j < (unsigned)mapiter->second.NumRows(); j++){
                    mapiter->second[j][i] = desiredCount.front();
                    desiredCount.pop();
                }
            }
        }
    }
    resultMatrix.kill();
    return;
}

void ecall_dispatcher::EncryptMat(vector<Ctxt>& ReturnedCtxts, vector<string>& Allleafmatrix)
{
    vector<Mat<size_t>> leafMatrix;
    size_t num_rawMat = Allleafmatrix.size();
    size_t num_subMat = num_rawMat / dim;
    size_t remainCols = num_rawMat % dim;
    Mat<size_t> Mattemp;
    cout << "Enclave: Check the Allleafmatrix size: " << num_rawMat << endl;
    if (num_rawMat < dim)
    {
        Mattemp.SetDims(dim, dim);
        MatInit(Mattemp);
        for (size_t j = 0; j < num_rawMat; ++j)
        {
            string temp = Allleafmatrix[j];
            for (size_t i = 0; i < /*rawdim = 13*/e_HTparams->numBits; ++i)
            {
                Mattemp[i][j] = temp[i] - '0';
            }
        }
        leafMatrix.push_back(Mattemp);
    }else{
        for (size_t k = 0; k < num_subMat; ++k)
        {
            Mattemp.SetDims(dim, dim);
            MatInit(Mattemp);
            size_t Idx = k * dim;
            for (size_t j = 0; j < dim; ++j)
            {
                string temp = Allleafmatrix[Idx++];
                for (size_t i = 0; i < /*rawdim = 13*/e_HTparams->numBits; ++i)
                {
                    Mattemp[i][j] = temp[i] - '0';
                }
            }
            leafMatrix.push_back(Mattemp);
            Mattemp.kill();
        }
        if (remainCols != 0)
        {
            Mattemp.SetDims(dim, dim);
            MatInit(Mattemp);
            size_t Idx = num_subMat * dim;
            for (size_t j = 0; j < remainCols; ++j)
            {
                string temp = Allleafmatrix[Idx++];
                for (size_t i = 0; i < /*rawdim = 13*/e_HTparams->numBits; ++i)
                {
                    Mattemp[i][j] = temp[i] - '0';
                }
            }
            leafMatrix.push_back(Mattemp);
        }
    }
    // encrypt ctxts in ReturnedCtxts
    for (size_t k = 0; k < leafMatrix.size(); ++k)
    {
        vector<long> cmsg(dim * dim);
        for(long i = 0; i < dim; ++i)
        {
            for(long j = 0; j < dim; ++j)
            {
                // size_t to long
                cmsg[i * dim + j] = leafMatrix[k][i][j];
            }
        }
        Ctxt ctemp(*activePubKey);
        e_context->getEA().encrypt(ctemp, *activePubKey, cmsg);
        ReturnedCtxts.push_back(ctemp);
    }
    // number of enclave ctxts
    num_returnedMat = ReturnedCtxts.size();
    cout << "Enclave: Check the ReturnedCtxts size: " << num_returnedMat << endl;
    return;
}

int ecall_dispatcher::HT_Classify(uint8_t* EvaCtxt, size_t EvaCtxt_len, size_t num_EvaCtxt)
{
    int ret = 0;

    string instance = string(EvaCtxt, EvaCtxt + EvaCtxt_len);
    Mat<size_t> instanceMatrix;
    instanceMatrix.SetDims(num_EvaCtxt, dim);
    MatInit(instanceMatrix);

    for (size_t i = 0; i < num_EvaCtxt; i++){
        for (size_t j = 0; j < EvaCtxt_len; j++){
            instanceMatrix[i][j] = instance[j]  - '0';
        }
    }
    cout << "Enclave: print instanceMatrix: " << endl;
    cout << instanceMatrix << endl;


    list<HoeffdingTree*>::iterator iter;
    // Collect leaf nodes path
    vector<string> pathCollector;
    // pair<countIndice, MajorityClass()>
    vector<pair<size_t, size_t>> labelCheck;
    for (iter = leafnodes.begin(); iter != leafnodes.end(); iter++)
    {
        if ((*iter)->remFeatures.empty())
        {
            // cout << "<count, majorityClass>: " << "<" << ((*iter)->CountIndice() - 1) << ", " << (*iter)->MajorityClass() << ">" << endl;
            labelCheck.push_back(make_pair(((*iter)->CountIndice() - 1), (*iter)->MajorityClass()));
        }else{
            // cout << "<count, majorityClass>: " << "<" << ((*iter)->CountIndice() - 2) << ", " << (*iter)->MajorityClass() << ">" << endl;
            labelCheck.push_back(make_pair(((*iter)->CountIndice() - 2), (*iter)->MajorityClass()));
        }
        pathCollector.push_back((*iter)->CurrentTreePath());
    }
    size_t num_leaf = pathCollector.size();

    // cout << "Enclave: Check TreeMat row size: " << e_HTparams->numBits << endl;
    // cout << "Enclave: Check TreeMat col size: " << num_leaf << endl;
    TreeMat.SetDims(dim, num_leaf);
    MatInit(TreeMat);
    for (size_t j = 0; j < num_leaf; j++) {
            string temp = pathCollector[j];
            for (size_t i = 0; i < e_HTparams->numBits; i++){
                TreeMat[i][j] = temp[i]  - '0';
            }
    }
    cout << "Enclave: print TreeMat: " << endl;
    cout << TreeMat << endl;

    // Plaintext Multiplication 
    Mat<size_t> predictMatrix;
    predictMatrix.SetDims(num_EvaCtxt, num_leaf);
    size_t multemp, addtemp;
    for (size_t i = 0; i < num_EvaCtxt; i++) {
        for (size_t j = 0; j < num_leaf; j++) {
            addtemp = 0;
            for (size_t k = 0; k < dim; k++) {
                multemp = instanceMatrix[i][k] * TreeMat[k][j];
                addtemp += multemp;
            }
            predictMatrix[i][j] = addtemp;
        }
    }
    cout << "Enclave: print predictMatrix: " << endl;
    cout << predictMatrix << endl;

    vector<size_t> predictions;
    for (size_t i = 0; i < num_EvaCtxt; i++) {
        for (size_t j = 0; j < num_leaf; j++) {
            if (predictMatrix[i][j] == labelCheck[j].first) {
                predictions.push_back(labelCheck[j].second);
                continue;
            }
        }
    }
    cout << "Enclave: print predictions: " << endl;
    cout << predictions << endl;

    // pathCollector.clear();
    // vector<string>().swap(pathCollector);

    // invoke the private function classify() for classification
    return ret;
}

void ecall_dispatcher::copy_arr_to_enclave(uint8_t* dst[], size_t num, uint8_t* src[], size_t lengths[]) {
  for (int i = 0; i < num; i++) {
    //size_t nlen = lengths[i];
    //check_host_buffer(src[i], nlen);
    dst[i] = (uint8_t*) malloc(lengths[i] * sizeof(uint8_t));
    memcpy(dst[i], src[i], lengths[i]);
  }
}

void ecall_dispatcher::free_array(uint8_t* arr[], size_t len) {
  for (int i = 0; i < len; i++) {
    free(arr[i]);
  }
}

void ecall_dispatcher::MatInit(Mat<size_t>& matInit)
{
    for (size_t i = 0; i < (unsigned)matInit.NumCols(); ++i){
        for (size_t j = 0; j < (unsigned)matInit.NumRows(); ++j){
            matInit[j][i] = 0;
        }
    }
}

void ecall_dispatcher::close()
{
    TRACE_ENCLAVE("Enclave: release context and HTparams!");
    // release context and keys
    delete e_context;
    delete e_HTparams;
    // activePubKey.reset();
    // activeSecKey.reset();
    TRACE_ENCLAVE("ecall_dispatcher::close");
}