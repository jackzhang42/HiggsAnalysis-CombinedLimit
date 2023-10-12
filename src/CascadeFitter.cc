#include "../interface/CascadeFitter.h"
// #include "../interface/RooMinimizerExt.h"
#include "../interface/ProfilingTools.h"
#include "../interface/utils.h"

#include <Math/MinimizerOptions.h>
#include <Math/IOptions.h>
#include <RooCategory.h>
#include <RooNumIntConfig.h>
#include <TStopwatch.h>
#include <RooStats/RooStatsUtils.h>

#include <iomanip>

boost::program_options::options_description CascadeFitter::options_("Cascade Fitter options");
std::vector<CascadeFitter::Algo> CascadeFitter::fallbacks_;
bool CascadeFitter::preScan_;
double CascadeFitter::approxPreFitTolerance_ = 0;
int CascadeFitter::approxPreFitStrategy_ = 0;
int  CascadeFitter::preFit_ = 0;
bool CascadeFitter::poiOnlyFit_;
bool CascadeFitter::singleNuisFit_;
bool CascadeFitter::setZeroPoint_ = true;
bool CascadeFitter::oldFallback_ = false;
bool CascadeFitter::firstHesse_ = false;
bool CascadeFitter::lastHesse_ = false;
int  CascadeFitter::minuit2StorageLevel_ = 0;
bool CascadeFitter::runShortCombinations = true;
float CascadeFitter::nuisancePruningThreshold_ = 0;
double CascadeFitter::discreteMinTol_ = 0.001;
std::string CascadeFitter::defaultMinimizerType_="Minuit2"; // default to minuit2 (not always the default !?)
std::string CascadeFitter::defaultMinimizerAlgo_="Migrad";
double CascadeFitter::defaultMinimizerTolerance_=1e-1;  
double CascadeFitter::defaultMinimizerPrecision_=-1.0;
int  CascadeFitter::strategy_=1; 

std::map<std::string,std::vector<std::string> > const CascadeFitter::minimizerAlgoMap_{
 {"Minuit"	 ,{"Migrad","Simplex","Combined","Scan"}}
,{"Minuit2" 	 ,{"Migrad","Simplex","Combined","Scan"}}
,{"GSLMultiMin"  ,{"ConjugateFR", "ConjugatePR", "BFGS", "BFGS2", "SteepestDescent"}}
};

CascadeFitter::CascadeFitter(RooAbsReal &nll, Mode mode, RooRealVar *poi) :
    nll_(nll),
    mode_(mode),
    //strategy_(0),
    poi_(poi),
    nuisances_(0),
    autoBounds_(false),
    poisForAutoBounds_(0),
    poisForAutoMax_(0)
{
    remakeMinimizer();
}

void CascadeFitter::remakeMinimizer() {
  std::cout << "this is a placeholder for CascadeFitter::remakeMinimizer()\n";
  // omitted CachingSimNLL things
  std::cout << "print nll details:\n";
  nll_.Print("v");
  minimizer_.reset();
  minimizer_.reset(new RooMinimizer(nll_));
}

bool CascadeFitter::freezeDiscParams(const bool freeze)
{
  std::cout << "this is a placeholder for CascadeFitter::remakeMinimizer()\n";
  return true;
}

void CascadeFitter::setAutoBounds(const RooArgSet *pois) 
{
    poisForAutoBounds_ = pois;
    autoBounds_ = (poisForAutoBounds_ != 0 || poisForAutoMax_ != 0);
}

void CascadeFitter::setAutoMax(const RooArgSet *pois) 
{
    poisForAutoMax_ = pois;
    autoBounds_ = (poisForAutoBounds_ != 0 || poisForAutoMax_ != 0);
}


bool CascadeFitter::improve(int verbose, bool cascade, bool forceResetMinimizer) 
{
  std::cout << "this is a placeholder for cascademinimizer::improve\n";
    if (forceResetMinimizer || !minimizer_.get()) remakeMinimizer();
    minimizer_->setPrintLevel(verbose-1);
   
    strategy_ = ROOT::Math::MinimizerOptions::DefaultStrategy(); // re-configure 

    minimizer_->setStrategy(strategy_);
    std::string nominalType(ROOT::Math::MinimizerOptions::DefaultMinimizerType());
    std::string nominalAlgo(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
    float       nominalTol(ROOT::Math::MinimizerOptions::DefaultTolerance());
    minimizer_->setEps(nominalTol);
    if (approxPreFitTolerance_ > 0) {
      double tol = std::max(approxPreFitTolerance_, 10. * nominalTol);
      do {
        if (verbose > 1) std::cout << "Running pre-fit with " << nominalType << "," << nominalAlgo << " and tolerance " << tol << std::endl;
        // Significance::MinimizerSentry minimizerConfig(nominalType+","+nominalAlgo, tol);
        minimizer_->setEps(tol);
        minimizer_->setStrategy(approxPreFitStrategy_);
        improveOnce(verbose-1, true);
        if (runtimedef::get("DBG_QUICKEXIT")) {
          exit(0);
        }
        minimizer_->setEps(nominalTol);
        minimizer_->setStrategy(strategy_);
      } while (autoBounds_ && !autoBoundsOk(verbose-1));
    }
    bool outcome;
    do {
      std::cout << "start improveonce\n";
      outcome = improveOnce(verbose-1);
      if (cascade && !outcome && !fallbacks_.empty()) {
        std::cout << "line 114 marker\n";
        int         nominalStrat(strategy_);
        if (verbose > 0) {
		std::cerr << "Failed minimization with " << nominalType << "," << nominalAlgo << " and tolerance " << nominalTol << std::endl;
		// Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Failed minimization with %s, %s and tolerance %g",__LINE__,nominalType.c_str(),nominalAlgo.c_str(),nominalTol)),Logger::kLogLevelDebug,__func__);
	}
        for (std::vector<Algo>::const_iterator it = fallbacks_.begin(), ed = fallbacks_.end(); it != ed; ++it) {
            // Significance::MinimizerSentry minimizerConfig(it->type + "," + it->algo, it->tolerance != Algo::default_tolerance() ? it->tolerance : nominalTol); // set the global defaults
            int myStrategy = it->strategy; if (myStrategy == Algo::default_strategy()) myStrategy = nominalStrat;
            if (nominalType != ROOT::Math::MinimizerOptions::DefaultMinimizerType() ||
                nominalAlgo != ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo() ||
                nominalTol  != ROOT::Math::MinimizerOptions::DefaultTolerance()     ||
                myStrategy  != nominalStrat) {
                if (verbose > 0) { 
			std::cerr << "Will fallback to minimization using " << it->algo << ", strategy " << myStrategy << " and tolerance " << it->tolerance << std::endl;
			// Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Will fallback to minimization using %s, strategy %d and tolerance %g",__LINE__,(it->algo).c_str(),myStrategy,it->tolerance)),Logger::kLogLevelDebug,__func__);
		}
                minimizer_->setEps(ROOT::Math::MinimizerOptions::DefaultTolerance());
                // minimizer_->setStrategy(myStrategy); 
                minimizer_->setStrategy(0);  // hardcoded
                outcome = improveOnce(verbose-2);
                if (outcome) break;
            }
        }
	
      }
    } while (autoBounds_ && !autoBoundsOk(verbose-1));


  return outcome;
}

bool CascadeFitter::improveOnce(int verbose, bool noHesse) 
{
  std::cout << "this is a placeholder for cascademinimizer::improveOnce\n";
    static int optConst = runtimedef::get("MINIMIZER_optimizeConst");
    static int rooFitOffset = runtimedef::get("MINIMIZER_rooFitOffset");
    std::string myType(ROOT::Math::MinimizerOptions::DefaultMinimizerType());
    std::string myAlgo(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
    int myStrategy = ROOT::Math::MinimizerOptions::DefaultStrategy();
    bool outcome = false;
    double tol = ROOT::Math::MinimizerOptions::DefaultTolerance();
    static int maxcalls = runtimedef::get("MINIMIZER_MaxCalls");
    if (!minimizer_.get()) remakeMinimizer();

    // freeze non active parameters if MINIMIZER_freezeDisassociatedParams enabled
    freezeDiscParams(true);

        // if (verbose+2>0) Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Minimisation configured with Type=%s, Algo=%s, strategy=%d, tolerance=%g",__LINE__,myType.c_str(),myAlgo.c_str(),myStrategy,tol)),Logger::kLogLevelInfo,__func__);
        // cacheutils::CachingSimNLL *simnll = setZeroPoint_ ? dynamic_cast<cacheutils::CachingSimNLL *>(&nll_) : 0;
        /* if (simnll) simnll->setZeroPoint();
         * if ((!simnll) && optConst) minimizer_->optimizeConst(std::max(0,optConst));
         * if ((!simnll) && rooFitOffset) minimizer_->setOffsetting(std::max(0,rooFitOffset)); */
        if (firstHesse_ && !noHesse) {
            minimizer_->setPrintLevel(std::max(0,verbose-3)); 
            minimizer_->hesse();
            // if (simnll) simnll->updateZeroPoint(); 
            minimizer_->setPrintLevel(verbose-1); 
        }
        int status = minimizer_->minimize("Minuit2", "Migrad");
        // int status = minimizer_->minimize(myType.c_str(), myAlgo.c_str());
        if (lastHesse_ && !noHesse) {
            // if (simnll) simnll->updateZeroPoint(); 
            minimizer_->setPrintLevel(std::max(0,verbose-3)); 
            status = minimizer_->hesse();
            minimizer_->setPrintLevel(verbose-1); 
          // if (verbose+2>0 ) Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Hesse finished with status=%d",__LINE__,status)),Logger::kLogLevelDebug,__func__);
        }
        // if (simnll) simnll->clearZeroPoint();
        outcome = (status == 0 || status == 1);
	if (status==1) std::cerr << "[WARNING] Minimisation finished with status 1 (covariance forced positive definite), this could indicate a problem with the minimim!" << std::endl;
/*     	if (verbose+2>0 ) {
 *     Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Minimisation finished with status=%d",__LINE__,status)),Logger::kLogLevelInfo,__func__);
 *     if (status==1) Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- finished with status 1 (covariance forced positive definite), this could indicate a problem with the minimim.",__LINE__)),Logger::kLogLevelDebug,__func__);
 *   }
 *     if (verbose+2>0 ){
 *      if  (outcome) Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Minimization success! status=0",__LINE__)),Logger::kLogLevelInfo,__func__);
 *      else Logger::instance().log(std::string(Form("CascadeMinimizer.cc: %d -- Minimization ended with latest status != 0 or 1",__LINE__)),Logger::kLogLevelDebug,__func__);
 *     }
 *  */
    // restore original params
    freezeDiscParams(false);


  return outcome;
}


bool CascadeFitter::minos(const RooArgSet & params , int verbose ) {
  std::cout << "this is a placeholder for cascademinimizer::minos\n";
  bool outcome = true;
  return outcome;
}

bool CascadeFitter::hesse(int verbose ) {
  std::cout << "this is a placeholder for cascademinimizer::hesse\n";
  bool outcome = true;
  return outcome;
}

bool CascadeFitter::iterativeMinimize(double &minimumNLL,int verbose, bool cascade){
  std::cout << "this is a placeholder for cascademinimizer::iterativeMinimize\n";
  bool outcome = true;
  return outcome;
}

bool CascadeFitter::minimize(int verbose, bool cascade) 
{
  std::cout << "this is a placeholder for cascademinimizer::minimize\n";
  std::cout << "currently testing index changing effect\n";
  // RooMinimizerExt minimizer(nll_);
  // minimizer.minimize("Minuit2", "Migrad");
  
  minimizer_->setPrintLevel(verbose-2);
  minimizer_->setStrategy(strategy_);

  RooArgSet reallyCleanParameters;
  std::unique_ptr<RooArgSet> nllParams(nll_.getParameters((const RooArgSet*)0));
  nllParams->remove(CascadeFitterGlobalConfigs::O().pdfCategories);
  // debug
  std::cout << "printing nllParams:\n";
  nllParams->Print("v");
  (nllParams)->snapshot(reallyCleanParameters); // should remove also the nuisance parameters from here!

  double minimumNLL = 10+nll_.getVal();
  std::vector<std::vector<bool>> contIndex;

  bool ret = true;
  
  multipleMinimize(reallyCleanParameters,ret,minimumNLL,verbose,cascade,0,contIndex);
  return ret;
}

bool CascadeFitter::multipleMinimize(const RooArgSet &reallyCleanParameters, bool& ret, double& minimumNLL, int verbose, bool cascade,int mode, std::vector<std::vector<bool> >&contributingIndeces){
    std::cout << "this is a placeholder for cascademinimizer::multipleMinimize\n";
    static bool freezeDisassParams = runtimedef::get(std::string("MINIMIZER_freezeDisassociatedParams"));
    static bool hideConstants = freezeDisassParams && runtimedef::get(std::string("MINIMIZER_multiMin_hideConstants"));
    static bool maskConstraints = freezeDisassParams && runtimedef::get(std::string("MINIMIZER_multiMin_maskConstraints"));
    static int maskChannels = freezeDisassParams ? runtimedef::get(std::string("MINIMIZER_multiMin_maskChannels")) : 0;
    // cacheutils::CachingSimNLL *simnll = dynamic_cast<cacheutils::CachingSimNLL *>(&nll_);

    //RooTrace::active(true);
    /* Different modes for minimization 
     Mode 0 -- Generate all combinations but only scan per-index 
     Mode 1 -- Generate only combinations which are orthogonal from the best fit after mode 0
	       Remove functions which cause increase in NLL > 10 (except best fit ones from previous mode)
     Mode 2 -- Full scan over the remaining combinations after mode 1
    */

    //std::cout << " At the start of the looping over the Indeces, minimum NLL is " << minimumNLL << std::endl; 
    // If the barlow-beeston minimisation is being used we can disable it temporarily,
    // saves time if we don't have to call enable/disable on the CMSHistErrorPropagators
    // repeatedly for no purpose
    int currentNoBarlowBeeston = runtimedef::get(std::string("MINIMIZER_no_analytic"));
    runtimedef::set("MINIMIZER_no_analytic", 1);
    
    double backupStrategy = ROOT::Math::MinimizerOptions::DefaultStrategy();
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);

    bool newDiscreteMinimum = false;

    RooArgList pdfCategoryIndeces = CascadeFitterGlobalConfigs::O().pdfCategories; 
    int numIndeces = pdfCategoryIndeces.getSize();
    
    // create all combinations of indeces 
    std::vector<int> pdfSizes;

    RooCategory *fPdf;

    std::vector<int> bestIndeces(numIndeces,0);

    // Set to the current best indeces
    for (int id=0;id<numIndeces;id++) {
	int c =((RooCategory*)(pdfCategoryIndeces.at(id)))->getIndex();
	bestIndeces[id]=c;
    } 

    if (mode==0) { // mode 0 makes the indeces
      contributingIndeces.clear();
      for (int id=0;id<numIndeces;id++){
    	int npdf = ((RooCategory*)(pdfCategoryIndeces.at(id)))->numTypes();
	std::vector<bool> indexFlags(npdf,true);
	contributingIndeces.push_back(indexFlags);
      }
    }

    // now find the number of available pdfs
    for (int id=0;id<numIndeces;id++){
    	int npdf = ((RooCategory*)(pdfCategoryIndeces.at(id)))->numTypes();
	pdfSizes.push_back(npdf);
    }

    // keep hold of best fitted parameters! 
    std::unique_ptr<RooArgSet> params;
    params.reset(nll_.getParameters((const RooArgSet *)0) );
    params->remove(CascadeFitterGlobalConfigs::O().pdfCategories);

    //take a snapshot of those parameters
    RooArgSet snap;
    params->snapshot(snap);
/* 
 *     if (maskChannels && simnll) {
 *         simnll->setMaskNonDiscreteChannels(true);
 *     }
 *     if (hideConstants && simnll) {
 *         simnll->setHideConstants(true);
 *         if (maskConstraints) simnll->setMaskConstraints(true);
 *         minimizer_.reset(); // will be recreated when needed by whoever needs it
 *     }
 *  */
    std::vector<std::vector<int> > myCombos;

    // Get All Permutations of pdfs
    if ( ( mode==0 ) /*&& runShortCombinations )*/ || mode ==1 ) myCombos = utils::generateOrthogonalCombinations(pdfSizes);
    else myCombos = utils::generateCombinations(pdfSizes);

    // Reorder to start from the "best indeces"
    //if (mode!=0) utils::reorderCombinations(myCombos,pdfSizes,bestIndeces);
    utils::reorderCombinations(myCombos,pdfSizes,bestIndeces);

    int numberOfCombinations = 1;
    if (mode==1 || mode==0) numberOfCombinations=myCombos.size();

    else {
      for (int i=0;i<numIndeces;i++){
   int nokpdfs=0;
         for (int j=0;j<pdfSizes[i];j++){
     nokpdfs+=contributingIndeces[i][j];
         }
   numberOfCombinations*=nokpdfs;
  }
    }

    std::vector<std::vector<int> >::iterator my_it = myCombos.begin();
    if (mode!=0) my_it++; // already did the best fit case
  
    TStopwatch tw; tw.Start();

    int fitCounter = 0;
    for (;my_it!=myCombos.end(); my_it++){

       bool isValidCombo = true;
  
       // int pdfIndex=0, changedIndex = -1;
       int pdfIndex = 0;
       int changedIndex = -1; 
       // Set the current indeces;
       std::vector<int> cit = *my_it;
       for (std::vector<int>::iterator it = cit.begin();
           it!=cit.end(); it++){

     isValidCombo &= (contributingIndeces)[pdfIndex][*it];
     if (!isValidCombo ) /*&& runShortCombinations)*/ continue;

          fPdf = (RooCategory*) pdfCategoryIndeces.at(pdfIndex);
          std::cout << "just to use changedIndex: " << changedIndex << std::endl;
          if (fPdf->getIndex() != *it) changedIndex = pdfIndex;
     fPdf->setIndex(*it);
     pdfIndex++;
       }
  
      if (!isValidCombo )/*&& runShortCombinations)*/ continue;
      
      if (verbose>2) {
  std::cout << "Setting indices := ";
  for (int id=0;id<numIndeces;id++) {
    std::cout << ((RooCategory*)(pdfCategoryIndeces.at(id)))->getIndex() << " ";
  }
        std::cout << std::endl;
      }

      if (fitCounter>0) params->assignValueOnly(reallyCleanParameters); // no need to reset from 0'th fit
/* 
 *       if (maskChannels == 2 && simnll) {
 *         for (int id=0;id<numIndeces;id++)  ((RooCategory*)(pdfCategoryIndeces.at(id)))->setConstant(id != changedIndex && changedIndex != -1);
 *         simnll->setMaskNonDiscreteChannels(true);
 *       } */
      // Remove parameters which are not associated to the current PDF (only works if using --X-rtd MINIMIZER_freezeDisassociatedParams)
      // freezeDiscParams(true);
/* 
 *       // FIXME can be made smarter than this
 *       if (mode_ == Unconstrained && poiOnlyFit_) {
 *         trivialMinimize(nll_, *poi_, 200);
 *       } */

      ret =  improve(verbose, cascade, freezeDisassParams);
/* 
 *       if (maskChannels == 2 && simnll) {
 *         for (int id=0;id<numIndeces;id++)  ((RooCategory*)(pdfCategoryIndeces.at(id)))->setConstant(false);
 *         simnll->setMaskNonDiscreteChannels(false);
 *       } */
      freezeDiscParams(false);


      fitCounter++;
      double thisNllValue = nll_.getVal();
      
      if ( thisNllValue < minimumNLL ){
    // Now we insert the correction ! 
                if (verbose>2) {
                    std::cout << " .... Found a better fit: new NLL = " << thisNllValue << " (improvement: " << (thisNllValue-minimumNLL) << std::endl;
                }
          minimumNLL = thisNllValue;	
                //std::cout << " .... Found a better fit! hoorah! " << minimumNLL << std::endl; 
        snap.assignValueOnly(*params);
    // set the best indeces again
    for (int id=0;id<numIndeces;id++) {
      if (bestIndeces[id] != ((RooCategory*)(pdfCategoryIndeces.at(id)))->getIndex() ) newDiscreteMinimum = true;
      bestIndeces[id]=((RooCategory*)(pdfCategoryIndeces.at(id)))->getIndex();	
    }
                if (verbose>2 && newDiscreteMinimum) {
                    std::cout << " .... Better fit corresponds to a new set of indices :=" ; 
                    for (int id=0;id<numIndeces;id++) { std::cout << " " << bestIndeces[id]; }
                    std::cout << std::endl;
                }
      }

      // FIXME this should be made configurable!
      double maxDeviation = 5;

      if (mode==1 )/*&& runShortCombinations)*/{

        if (thisNllValue > minimumNLL+maxDeviation){
    // Step 1, find out which index just changed 
    int modid   =0;
    int modcount=0;

          for (int id=0;id<numIndeces;id++) {
      RooCategory* thisCat = (RooCategory*)(pdfCategoryIndeces.at(id));
      if (thisCat->getIndex()!=bestIndeces[id]){
        modid=id;
        modcount++;
      }
    }
    
    if (modcount==1){
      // Step 2, remove its current index from the allowed indexes
      RooCategory* thisCat = (RooCategory*)(pdfCategoryIndeces.at(modid));
      int cIndex = thisCat->getIndex();
      if (cIndex!=bestIndeces[modid]){ // don't remove the best pdf for this index!
      (contributingIndeces)[modid][cIndex]=false;
      }
    }
        }
     }
 
    }
/* 
 *     // Assign best values ;
 *     for (int id=0;id<numIndeces;id++) {
 *   ((RooCategory*)(pdfCategoryIndeces.at(id)))->setIndex(bestIndeces[id]);	
 *     } 
 *     params->assignValueOnly(snap);
 * 
 *     runtimedef::set("MINIMIZER_no_analytic", currentNoBarlowBeeston);
 *     ROOT::Math::MinimizerOptions::SetDefaultStrategy(backupStrategy);
 * 
 *     tw.Stop(); if (verbose > 2) std::cout << "Done " << myCombos.size() << " combinations in " << tw.RealTime() << " s. New discrete minimum? " << newDiscreteMinimum << std::endl;
 * 
 *     if (maskChannels && simnll) {
 *         simnll->setMaskNonDiscreteChannels(false);
 *     }
 *     if (hideConstants && simnll) {
 *         simnll->setHideConstants(false);
 *         if (maskConstraints) simnll->setMaskConstraints(false);
 *         minimizer_.reset(); // will be recreated when needed by whoever needs it
 *     }
 *   */
 
    return newDiscreteMinimum;
}

void CascadeFitter::initOptions() 
{
  std::cout << "this is a placeholder for cascademinimizer::initOptions\n";
}

bool CascadeFitter::checkAlgoInType(std::string type, std::string algo){
  std::cout << "this is a placeholder for cascademinimizer::checkAlgoInType\n";
  bool outcome = true;
  return outcome;
}

void CascadeFitter::applyOptions(const boost::program_options::variables_map &vm) 
{
  std::cout << "this is a placeholder for cascademinimizer::applyOptions\n";
}

void CascadeFitter::trivialMinimize(const RooAbsReal &nll, RooRealVar &r, int points) const {
  std::cout << "this is a placeholder for cascademinimizer::trivialMinimize\n";
}

bool CascadeFitter::autoBoundsOk(int verbose) {
  std::cout << "this is a placeholder for cascademinimizer::autoBoundsOk\n";
  bool outcome = true;
  return outcome;
}
