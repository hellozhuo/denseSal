//04 2017

#include <iostream>
#include <vector>
#include<fstream>

#include <CmCode/CmLib/CmLib.h>
#include<opencv2/opencv.hpp>
#include "DenseCRF.h"
#include "mexopencv.hpp"

void getParaFromMat(std::vector<float>& paras, const cv::Mat& paraM)
{
    assert(paraM.type()==CV_64F);
    double* parabuf = (double*)paraM.data;
    const int N = paras.size();
    for(int i = 0; i<N; i++)
    {
        paras[i] = parabuf[i];
    }
}

void denseCut(const cv::Mat& fgProb, cv::Mat& res,
	const cv::Mat& normLab, const cv::Mat& srcRGB,
    const CmGMM& _fGMM, const CmGMM& _bGMM, int iter, double wei,
	float w1, float w2, float w3, float alpha, float beta, float gama, float mu)
{
	//6, 10, 2, 20, 33, 3, 41
	const int _h = srcRGB.rows;
	const int _w = srcRGB.cols;
	cv::Mat _unary2f(_h, _w, CV_32FC2);
	DenseCRF2D _crf(_w, _h, 2);
	res.create(_h, _w, CV_64F);

	if (w1 != 0)
		_crf.addPairwiseBilateral(alpha, alpha, beta, beta, beta, srcRGB.data, w1);
	if (w2 != 0)
		_crf.addPairwiseGaussian(gama, gama, w2);
	if (w3 != 0)
		_crf.addPairwiseColorGaussian(mu, mu, mu, srcRGB.data, w3);

	cv::Vec2f* unryV = (cv::Vec2f*)_unary2f.data;
	cv::Vec3f* labbuf = (cv::Vec3f*)normLab.data;
	float backW = 1.0;
    float weif = wei;
	for (int i = 0; i < _w*_h; i++)
	{
		float foreP = _fGMM.P(labbuf[i]), backP = _bGMM.P(labbuf[i]);
		float prb = weif * foreP / (foreP + backW*backP + 1e-8f);
// 		unryV[i] = cv::Vec2f(prb, 1 - prb);
        unryV[i] = cv::Vec2f(-fast_log(1-prb+FLT_EPSILON),-fast_log(prb+FLT_EPSILON));
	}

	_crf.setUnaryEnergy(_unary2f.ptr<float>(0));
	float* prob = _crf.binarySeg((float*)fgProb.data, iter, 1.f);
	double* resbuf = (double*)res.data;
	const int N = _w * _h;
	for (int i = 0; i<N; i++, prob += 2)
		resbuf[i] = prob[1] / (prob[0] + prob[1] + 1e-20f);
}

void getFBprobFromGMM(cv::Mat& res,
	const cv::Mat& srcRGB, const cv::Mat& GMM_fg, const cv::Mat& GMM_bg, const cv::Mat& fgProb,
        int iter, double wei, const std::vector<float>& paras)
{
	assert(srcRGB.type() == CV_8UC3 && GMM_bg.type() == CV_32F && GMM_fg.type() == CV_32F);
    assert(fgProb.type() == CV_32F);
	assert(srcRGB.size() == GMM_bg.size() && srcRGB.size() == GMM_fg.size());

	cv::Mat normLab;
	CmGMM _bGMM(5), _fGMM(5); // Background and foreground GMM
	Mat _bGMMidx1i, _fGMMidx1i;	// Background and foreground GMM components, supply memory for GMM, not used for Grabcut 

	cv::cvtColor(srcRGB, normLab, CV_RGB2Lab);
	normLab.convertTo(normLab, CV_32F, 1 / 255.0);

	_bGMM.BuildGMMs(normLab, _bGMMidx1i, GMM_bg);
	_fGMM.BuildGMMs(normLab, _fGMMidx1i, GMM_fg);

	_bGMM.RefineGMMs(normLab, _bGMMidx1i, GMM_bg);
	_fGMM.RefineGMMs(normLab, _fGMMidx1i, GMM_fg);

	denseCut(fgProb, res, normLab, srcRGB, _fGMM, _bGMM, iter, wei, paras[0],paras[1],paras[2],paras[3],paras[4],paras[5],paras[6]);
    
}

/**
* Main entry called from Matlab
* @param nlhs number of left-hand-side arguments
* @param plhs pointers to mxArrays in the left-hand-side
* @param nrhs number of right-hand-side arguments
* @param prhs pointers to mxArrays in the right-hand-side
*/

//dst = getFilteredLogFg(srcRGB, GMM_fg, GMM_bg, fgProb, iter=2, wei = 0.8, paras = [...]);
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    // Check the number of arguments
    if (nrhs < 4 || nrhs > 7 || nlhs != 1)
        mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
    
    // Decide second arguments
    int iter = 2;
    double wei = 0.8;
    std::vector<float> paras = {6.0, 10.0, 2.0, 20.0, 33.0, 3.0, 41.0};
    if(nrhs >= 5)
        iter = ((MxArray)prhs[4]).toInt();
    if(nrhs >= 6)
        wei = ((MxArray)prhs[5]).toDouble();
    if(nrhs >= 7)
    {
        cv::Mat paraM = ((MxArray)prhs[6]).toMat();
        getParaFromMat(paras, paraM);
    }

	cv::Mat srcRGB(((MxArray)prhs[0]).toMat());
	cv::Mat GMM_fg(((MxArray)prhs[1]).toMat());
	cv::Mat GMM_bg(((MxArray)prhs[2]).toMat());
    cv::Mat fgProb(((MxArray)prhs[3]).toMat());

	cv::Mat dst;
    
    // Apply
	getFBprobFromGMM(dst, srcRGB, GMM_fg, GMM_bg, fgProb, iter, wei,paras);
    plhs[0] = MxArray(dst);
}
