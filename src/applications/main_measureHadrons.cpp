/*
 * main_measureHadrons.cpp
 *
 * Luis Altenkort, 6 Jan 2021
 *
 * Main application for measuring hadronic correlators.
 *
 */

#include "../simulateqcd.h"
#include "../modules/measureHadrons/measureHadrons.h"

#define USE_GPU true
#if SINGLEPREC
#define PREC float
#else
#define PREC double
#endif

template< class GaugeField >
void readConfFormat( GaugeField & gauge , measureHadronsParam<PREC>& lp ){
    if ( lp.format() == "nersc" ){
        gauge.readconf_nersc( lp.GaugefileName() );
    }
    else if ( lp.format() == "milc" ){
        gauge.readconf_milc( lp.GaugefileName() );
    }
    else if ( lp.format() == "ildg" ){
        gauge.readconf_ildg( lp.GaugefileName() );
    }
    else {
        throw std::runtime_error(stdLogger.fatal("Unknown conf format use one of {nersc , milc , ildg }."));
    }

}
int main(int argc, char *argv[]) {
    try {
        stdLogger.setVerbosity(TRACE);

        CommunicationBase commBase(&argc, &argv);
        measureHadronsParam<PREC> lp;
        lp.readfile(commBase, "../parameter/tests/measureHadronsTest.param", argc, argv);
        commBase.init(lp.nodeDim());
        lp.check_for_nonsense();

        const size_t HaloDepth = 2; //! reason: HISQ Dslash
        const size_t HaloDepthSpin = 4; //! reason: HISQ Dslash
        const size_t NStacks = 1; //TODO add support for multiple sources at the same time i.e. nstacks>1
        const CompressionType compHISQ = R18 ;
        const CompressionType compstdStag = R14 ;
        
        initIndexer(HaloDepth, lp, commBase);
        
        if ( lp.action() == "HISQ" )
        {
            Gaugefield<PREC, USE_GPU, HaloDepth , compHISQ> gauge(commBase);
            if (lp.use_unit_conf()){
                rootLogger.info("Using unit configuration for tests/benchmarks");
                gauge.one();
            } else {
                rootLogger.info("Read configuration");
                //gauge.readconf_nersc(lp.GaugefileName());
                readConfFormat< Gaugefield<PREC, USE_GPU, HaloDepth , compHISQ> >( gauge , lp );
            }
            gauge.updateAll();
    
            //! Check plaquette
            GaugeAction<PREC, USE_GPU, HaloDepth , compHISQ > gAction(gauge);
            PREC plaq;
            plaq = gAction.plaquette();
            rootLogger.info( "plaquette: " ,  plaq ) ;

        
            measureHadrons<PREC, USE_GPU, HaloDepth, HaloDepthSpin, Even, NStacks, compHISQ, compstdStag> mesons(commBase, lp );
            mesons.compute_HISQ_correlators( gauge) ;
            mesons.write_correlators_to_file();
        }
        else if ( lp.action() == "stdStag" )
        {
            Gaugefield<PREC, USE_GPU, HaloDepth , compstdStag> gauge(commBase);
            if (lp.use_unit_conf()){
                rootLogger.info( "Using unit configuration for tests/benchmarks" );
                gauge.one();
            } else {
                rootLogger.info( "Read configuration" );
                //gauge.readconf_nersc(lp.GaugefileName());
                readConfFormat< Gaugefield<PREC, USE_GPU, HaloDepth , compstdStag> >( gauge , lp );

            }
            gauge.updateAll();

            //! Check plaquette
            GaugeAction<PREC, USE_GPU, HaloDepth , compstdStag > gAction(gauge);
            PREC plaq;
            plaq = gAction.plaquette();
            rootLogger.info("plaquette: " ,  plaq);

        
            measureHadrons<PREC, USE_GPU, HaloDepth, HaloDepthSpin, Even, NStacks, compHISQ, compstdStag> mesons(commBase, lp );
            mesons.compute_stdStag_correlators( gauge) ;
            mesons.write_correlators_to_file();
        }
        else{
            throw std::runtime_error(stdLogger.fatal("Unknown action! Choose one from {HISQ , stdStag}."));
        }

        //!Additional features we want to have:
        //TODO create a test routine that checks if the program works for all possible parameters on a small configuration
        //     using results from the BielefeldGPUCode as a reference
        //TODO add the option for multiple sources (maybe first one after another and then in parallel using stacked
        // spinorfields and multiRHS inverter?)
        //TODO add an option to choose the boundary conditions for the fermions (?)
    }

    catch (const std::runtime_error &error) {
        return -1;
    }
    return 0;
}
