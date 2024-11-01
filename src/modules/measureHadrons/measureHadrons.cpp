//
// Created by luis on 11.01.21.
//

#include "measureHadrons.h"

//! Functor to compute contractions for stacked spinor with nstacks=3
template<class floatT, size_t HaloDepth>
struct contractPropagators {

    Vect3arrayAcc<floatT> _acc_i;
    Vect3arrayAcc<floatT> _acc_j;

    contractPropagators(Spinorfield<floatT, false, All, HaloDepth, 3> &spinor_i,
                        Spinorfield<floatT, false, All, HaloDepth, 3> &spinor_j) :
            _acc_i(spinor_i.getAccessor()), _acc_j(spinor_j.getAccessor()) {}


    __host__ __device__ GPUcomplex<floatT> operator()(gSite site) {
        gSiteStack red = GIndexer<All, HaloDepth>::getSiteStack(site, 0);
        gSiteStack green = GIndexer<All, HaloDepth>::getSiteStack(site, 1);
        gSiteStack blue = GIndexer<All, HaloDepth>::getSiteStack(site, 2);
        GPUcomplex<floatT> contraction =
                  dot_prod<floatT>(_acc_i.template getElement<floatT>(red  ), _acc_j.template getElement<floatT>(red))
                + dot_prod<floatT>(_acc_i.template getElement<floatT>(green), _acc_j.template getElement<floatT>(green))
                + dot_prod<floatT>(_acc_i.template getElement<floatT>(blue ), _acc_j.template getElement<floatT>(blue));
        return contraction;
    }
};

template<class floatT, bool onDevice, size_t HaloDepth, size_t HaloDepthSpin, Layout Source,  size_t NStacks, CompressionType compHISQ , CompressionType compstdStag >
template<const size_t n_color>
void measureHadrons<floatT, onDevice, HaloDepth, HaloDepthSpin, Source, NStacks, compHISQ, compstdStag>::setup_memory_StaggeredType< n_color >(
    std::vector<Spinorfield<floatT,false,All,0,n_color>> & prop_containers
    ){
    for(size_t i = 0; i < _n_masses; i++) {
        prop_containers.emplace_back(std::move(Spinorfield<floatT,false,All,0,n_color>(_commBase, "Spinorfield"+_lp.mass_labels.get()[i])));
    }

};

template<class floatT, bool onDevice, size_t HaloDepth, size_t HaloDepthSpin, Layout Source,  size_t NStacks, CompressionType compHISQ , CompressionType compstdStag >
template<const size_t n_color>
void measureHadrons<floatT, onDevice, HaloDepth, HaloDepthSpin, Source, NStacks, compHISQ, compstdStag>::contract_prop_containers< n_color >(
    std::vector<Spinorfield<floatT,false,All,0,n_color>> & prop_containers
    ){
    //! contract the prop_containers for all mass combinations. look into the functor "contractPropagators" at the top
    //! to see what's going on
    size_t mass_index = 0;
    for (int m_i = 0; m_i < _n_masses; m_i++) {
        for (int m_j = 0; m_j <= m_i; m_j++) {
            contractPropagators<floatT, false> contractProp(prop_containers.at(m_i),prop_containers.at(m_j));
            _contracted_propagators.at(mass_index).template iterateOverBulk<All,0>(contractProp);
            mass_index++;
        }
    }
};

template< class floatT, bool onDevice, size_t HaloDepth, size_t HaloDepthSpin, Layout Source,  size_t NStacks, CompressionType compHISQ, CompressionType compstdStag > 
template< class StaggeredTypeAction , const size_t n_color >
void measureHadrons<floatT, onDevice, HaloDepth, HaloDepthSpin, Source, NStacks, compHISQ, compstdStag>::invert_StaggeredType_Dslash<StaggeredTypeAction , n_color>( 
    StaggeredTypeAction & dslash , 
    size_t massIndex , 
    floatT mass ,
    std::vector<Spinorfield<floatT,false,All,0,n_color>>  & prop_containers ){
    
    floatT cg_residue = _individual_cg_residues[massIndex];
    int cg_max_iter = _individual_cg_max_iters[massIndex];

    rootLogger.info( "cg_residue = " , cg_residue , ", max_iter = " , cg_max_iter );

    //! host
    Spinorfield<floatT, false, Even, HaloDepthSpin, NStacks> h_Ge(_commBase); //! even part of propagator
    Spinorfield<floatT, false, Odd, HaloDepthSpin, NStacks> h_Go(_commBase); //! odd part of propagator
    //! device
    Spinorfield<floatT, true, Even, HaloDepthSpin, NStacks> d_Se(_commBase); //! even source
    Spinorfield<floatT, true, Even, HaloDepthSpin, NStacks> d_Ge(_commBase); //! even part of propagator
    Spinorfield<floatT, true, Odd, HaloDepthSpin, NStacks> d_Go(_commBase); //! odd part of propagator

    ConjugateGradient<floatT,NStacks> cg;
    Vect3arrayAcc<floatT> h_Ge_acc = h_Ge.getAccessor();               //! host accessor for even propagator
    Vect3arrayAcc<floatT> h_Go_acc = h_Go.getAccessor();               //! host accessor for odd propagator
 
    //! first color
    const int color0 = 0;
    d_Se.setPointSource(_current_source, color0, mass);
 
    //! Solve Ax=b for x. (i.e. find A^{-1} and multiply it from the left)
    //! A: (m^2-D_eo*D_oe)      [sometimes sloppily called "MdaggM"]
    //! x: G_e                  [even part of propagator]
    //! b: m*S_e                [point source]
    cg.invert(dslash, d_Ge, d_Se, cg_max_iter, cg_residue); //! this takes up most of the computation time
    h_Ge = d_Ge;
    //! iterate over solution in h_Ge and copy its contents into the correct places in prop_container[i]
    prop_containers[massIndex].template iterateOverEvenBulkAtStack<color0>(
            returnSpinor<floatT, false, Even, HaloDepthSpin, NStacks>(h_Ge));
 
    //! Now we want to compute G_o = -1/m * D_oe G_e:
    //! First apply D_oe to G_e, then factor in -1/m
    dslash.Dslash(d_Go, d_Ge, true);
    d_Go = (static_cast<floatT>(-1.) / mass) * d_Go;
    h_Go = d_Go;
    //! same for odd part
    prop_containers[massIndex].template iterateOverOddBulkAtStack<color0>(
            returnSpinor<floatT, false, Odd, HaloDepthSpin, NStacks>(h_Go));
 
    //! second color
    const int color1 = 1;
    d_Se.setPointSource(_current_source, color1, mass);
    cg.invert(dslash, d_Ge, d_Se, cg_max_iter, cg_residue); //! this takes up most of the computation time
    h_Ge = d_Ge;
    prop_containers[massIndex].template iterateOverEvenBulkAtStack<color1>(
            returnSpinor<floatT, false, Even, HaloDepthSpin, NStacks>(h_Ge));
    dslash.Dslash(d_Go, d_Ge, true);
    d_Go = (static_cast<floatT>(-1.) / mass) * d_Go;
    h_Go = d_Go;
    prop_containers[massIndex].template iterateOverOddBulkAtStack<color1>(
            returnSpinor<floatT, false, Odd, HaloDepthSpin, NStacks>(h_Go));
 
    //! third color
    const int color2 = 2;
    d_Se.setPointSource(_current_source, color2, mass);
    cg.invert(dslash, d_Ge, d_Se, cg_max_iter, cg_residue); //! this takes up most of the computation time
    h_Ge = d_Ge;
    prop_containers[massIndex].template iterateOverEvenBulkAtStack<color2>(
            returnSpinor<floatT, false, Even, HaloDepthSpin, NStacks>(h_Ge));
    dslash.Dslash(d_Go, d_Ge, true);
    d_Go = (static_cast<floatT>(-1.) / mass) * d_Go;
    h_Go = d_Go;
    prop_containers[massIndex].template iterateOverOddBulkAtStack<color2>(
            returnSpinor<floatT, false, Odd, HaloDepthSpin, NStacks>(h_Go));


}

template<class floatT, bool onDevice, size_t HaloDepth, size_t HaloDepthSpin, Layout Source,  size_t NStacks, CompressionType compHISQ , CompressionType compstdStag >
void measureHadrons<floatT, onDevice, HaloDepth, HaloDepthSpin, Source, NStacks, compHISQ, compstdStag>::compute_StaggeredType_correlator(){
    size_t mass_index = 0;
    for (size_t corr_axis_index = 0; corr_axis_index < _n_correlator_axes ; corr_axis_index++ ) 
    {
    
    
    rootLogger.info("Computing correlators along " ,  _lp.correlator_axes.get()[corr_axis_index] ,  "-axis.");
        for (size_t momentum_index = 0; momentum_index < _momentum_elements_per_axis[corr_axis_index]; momentum_index++)
        {
            floatT p1 = static_cast<floatT>(2 * M_PI * _momenta[ _momenta_start_index[corr_axis_index] + 3 * momentum_index]     / _global_lat_extents[_axes_indices[ 4 * corr_axis_index + 1]]) ;
            floatT p2 = static_cast<floatT>(2 * M_PI * _momenta[ _momenta_start_index[corr_axis_index] + 3 * momentum_index + 1 ] / _global_lat_extents[_axes_indices[ 4 * corr_axis_index + 2]]) ;
            floatT p3 = static_cast<floatT>(2 * M_PI * _momenta[ _momenta_start_index[corr_axis_index] + 3 * momentum_index + 2 ] / _global_lat_extents[_axes_indices[ 4 * corr_axis_index + 3]]) ; 
            rootLogger.info(  _lat_extents[_axes_indices[ 4 * corr_axis_index + 1]] , " " , _lat_extents[_axes_indices[ 4 * corr_axis_index + 2]] , " " , _lat_extents[_axes_indices[ 4 * corr_axis_index + 3]] ) ;
            rootLogger.info(  _global_lat_extents[_axes_indices[ 4 * corr_axis_index + 1]] , " " , _global_lat_extents[_axes_indices[ 4 * corr_axis_index + 2]] , " " , _global_lat_extents[_axes_indices[ 4 * corr_axis_index + 3]] ) ;

            //rootLogger.info(  _momenta[_momenta_start_index[corr_axis_index] + 3 * momentum_index] ," ", _momenta[_momenta_start_index[corr_axis_index] + 3 * momentum_index + 1] ," ", _momenta[_momenta_start_index[corr_axis_index] + 3 * momentum_index + 2] );
            //rootLogger.info(  p1 ," ", p2 ," ", p3 );
            std::complex<floatT> ImaginaryUnit(0.0,1.0) ;      
        
        //! Reduce contracted propagators and project to different quantum numbers (channels).
        //! We do this here without functor syntax because then we can easily factor in the phase factors while reducing.
            for (int channel = 1; channel <= 8; channel++){
                rootLogger.info("Project to M" ,  channel);
                mass_index = 0;
                for (int m_i = 0; m_i < _n_masses; m_i++){
                    for (int m_j = 0; m_j <= m_i; m_j++){
                        LatticeContainerAccessor prop_acc = _contracted_propagators[mass_index].getAccessor();
        
                        for (int w = 0; w < _corr_ls[corr_axis_index] ; w++) {
                            GPUcomplex<floatT> result(0.0,0.0);
        
                            //! reduce over the volume without correlator axis
                            for (int r = 0; r < _lat_extents[_axes_indices[ 4 * corr_axis_index + 1]]; ++r) {
                                for (int u = 0; u < _lat_extents[_axes_indices[ 4 * corr_axis_index + 2 ]]; ++u) {
                                    for (int v = 0; v < _lat_extents[_axes_indices[ 4 * corr_axis_index + 3 ]]; ++v) {
                                        sitexyzt coord(0,0,0,0);
                                        coord[_axes_indices[4 * corr_axis_index + 0]] = w;
                                        coord[_axes_indices[4 * corr_axis_index + 1]] = r;
                                        coord[_axes_indices[4 * corr_axis_index + 2]] = u;
                                        coord[_axes_indices[4 * corr_axis_index + 3]] = v;
                                        sitexyzt globalCoord = GIndexer<All, HaloDepthSpin>::getLatData().globalPos(coord);
                                        
                                        std::complex<floatT> ValueExp = std::exp( ImaginaryUnit * static_cast<floatT>( 
                                            (p1 * globalCoord[_axes_indices[4 * corr_axis_index + 1]] 
                                           + p2 * globalCoord[_axes_indices[4 * corr_axis_index + 2]] 
                                           + p3 * globalCoord[_axes_indices[4 * corr_axis_index + 3]] ))) ;
                                        GPUcomplex<floatT> ValueExp2(ValueExp.real(),ValueExp.imag());
                                        
                                        
                                        int isite = GIndexer<All, HaloDepthSpin>::coordToIndex_Bulk(coord);
        
                                        //! it's okay that we're using the sublattice coordinates for the phase_factors, since the (sub)lattice extents are always even.
                                        result += prop_acc.template getElement<GPUcomplex<floatT>>(isite) 
                                                * static_cast<floatT>(phase_factor(channel, r, u, v))
                                                * ValueExp2 ;
                                    }
                                }
                            }
                            //! reduce the result from multiple gpus
                            _correlators[corr_index(( w - _current_source[_axes_indices[4 * corr_axis_index + 0]] + _corr_ls[corr_axis_index])% _corr_ls[corr_axis_index],
                                                    channel, 
                                                    mass_index , 
                                                    corr_axis_index ,
                                                    momentum_index  )] += _commBase.reduce(result);
                        }
                        mass_index++;
                    }
                }
            }
        }
    }
    rootLogger.info("Done calculating correlators");


}

template<class floatT, bool onDevice, size_t HaloDepth, size_t HaloDepthSpin, Layout Source,  size_t NStacks, CompressionType compHISQ , CompressionType compstdStag >
void measureHadrons<floatT, onDevice, HaloDepth, HaloDepthSpin, Source, NStacks, compHISQ, compstdStag>::compute_HISQ_correlators(Gaugefield<floatT,onDevice,HaloDepth,compHISQ>& gauge){

    //! setup memory for propagator results
    const size_t n_color = 3;
    std::vector<Spinorfield<floatT,false,All,0,n_color>> prop_containers;

    setup_memory_StaggeredType( prop_containers ) ;
    
    
    for(int i = 0; i < _n_masses; i++){
        floatT mass = _lp.masses.get()[i];
        
        
            floatT naik_epsilon = _individual_naik_epsilons[i];
    
            rootLogger.info( "------------------------------------------------------------------------" ) ;
            rootLogger.info( "Computing propagators: m_" , _lp.mass_labels.get()[i] , " = " , mass , ", naik_eps = " , naik_epsilon , ", " ) ;
    
    
            //! Apply HISQ smearing and perform staggered transformation
            Gaugefield<floatT, onDevice, HaloDepth, R18> gauge_smeared(_commBase);
            Gaugefield<floatT, onDevice, HaloDepth, U3R14> gauge_Naik(_commBase);
            HisqSmearing<floatT, onDevice, HaloDepth,compHISQ,R18,R18,U3R14> smearing(gauge, gauge_smeared, gauge_Naik, naik_epsilon);
            smearing.SmearAll( 0.0 , true );
    
            HisqDSlash<floatT, onDevice, Even, HaloDepth, HaloDepthSpin, NStacks> dslash(gauge_smeared, gauge_Naik,
                                                                                         mass, naik_epsilon);
    
            invert_StaggeredType_Dslash< HisqDSlash<floatT, onDevice, Even, HaloDepth, HaloDepthSpin, NStacks> , n_color>( 
            dslash , 
            i , 
            mass , 
            prop_containers ) ;
            
            
   
    }

    
    contract_prop_containers( prop_containers ) ;
    compute_StaggeredType_correlator() ;
    

}

template<class floatT, bool onDevice, size_t HaloDepth, size_t HaloDepthSpin, Layout Source,  size_t NStacks, CompressionType compHISQ , CompressionType compstdStag >
void measureHadrons<floatT, onDevice, HaloDepth, HaloDepthSpin, Source, NStacks, compHISQ, compstdStag>::compute_stdStag_correlators(Gaugefield<floatT,onDevice,HaloDepth,compstdStag>& gauge){

    //! setup memory for propagator results
    const size_t n_color = 3;
    std::vector<Spinorfield<floatT,false,All,0,n_color>> prop_containers;

    setup_memory_StaggeredType( prop_containers ) ;
    
    
    for(int i = 0; i < _n_masses; i++){
        floatT mass = _lp.masses.get()[i];
        


        rootLogger.info( "------------------------------------------------------------------------" );
        rootLogger.info( "Computing propagators: m_" , _lp.mass_labels.get()[i] , " = " , mass  , ", ") ;



        stdStagDSlash<floatT, onDevice, Even, HaloDepth, HaloDepthSpin, NStacks> dslash( gauge , mass );

        invert_StaggeredType_Dslash< stdStagDSlash<floatT, onDevice, Even, HaloDepth, HaloDepthSpin, NStacks> , n_color>( 
        dslash , 
        i , 
        mass , 
        prop_containers ) ;
        
  
    }

    
    contract_prop_containers( prop_containers ) ;
    compute_StaggeredType_correlator() ;
    
}

template<class floatT, bool onDevice, size_t HaloDepth, size_t HaloDepthSpin, Layout Source, size_t NStacks, CompressionType compHISQ, CompressionType compstdStag>
void measureHadrons<floatT, onDevice, HaloDepth, HaloDepthSpin, Source, NStacks, compHISQ, compstdStag>::write_correlators_to_file() {

    std::map<int, std::string> name_mapping = {
            {1, "scalar "},
            {2, "psdsclr"},
            {3, "axvctrx"},
            {4, "axvctry"},
            {5, "axvctrt"},
            {6, "vectorx"},
            {7, "vectory"},
            {8, "vectort"}
    };

    
    std::stringstream full_path_to_file;
    for (size_t corr_axis_index = 0; corr_axis_index < _n_correlator_axes ; corr_axis_index++ )
    {

        full_path_to_file.str("") ;
        std::array<std::string,4> corr_axis_prefixes = {"x", "y", "z", "t"};
        std::string corr_axis_prefix = corr_axis_prefixes[_axes_indices[4 * corr_axis_index ]];
        std::string momentum_prefix1 = corr_axis_prefixes[_axes_indices[4 * corr_axis_index + 1 ]];
        std::string momentum_prefix2 = corr_axis_prefixes[_axes_indices[4 * corr_axis_index + 2 ]];
        std::string momentum_prefix3 = corr_axis_prefixes[_axes_indices[4 * corr_axis_index + 3 ]];
        
        name_mapping[3] = "axvctr" + momentum_prefix1 ;
        name_mapping[6] = "vector" + momentum_prefix1 ;
        name_mapping[4] = "axvctr" + momentum_prefix2 ;
        name_mapping[7] = "vector" + momentum_prefix2 ;
        name_mapping[5] = "axvctr" + momentum_prefix3 ;
        name_mapping[8] = "vector" + momentum_prefix3 ;


        full_path_to_file << _lp.measurements_dir() << "/" << "corr_" << _lp.action() << "_" << _lp.source_type() << "_" << corr_axis_prefix << "_"
                          << _current_source.x << "-" << _current_source.y << "-" << _current_source.z << "-"
                          << _current_source.t <<  _lp.fileExt();
        FileWriter file(_commBase, _lp, full_path_to_file.str());
        rootLogger.info("Writing correlators to file: " ,  full_path_to_file.str());
        LineFormatter header = file.header();
        header << corr_axis_prefix+"  k"+momentum_prefix1+"  k"+momentum_prefix2+"  k"+momentum_prefix3 << "real " << "imag ";
        header.endLine();

        #define format std::fixed << std::setprecision(6)

        size_t mass_index = 0;
        for (int m_i = 0; m_i < _n_masses; m_i++){
            for (int m_j = 0; m_j <= m_i; m_j++){
                file << "# m_" << _lp.mass_labels[m_i] << "=" << format << _lp.masses[m_i]
                     << ", eps=" << format << _individual_naik_epsilons[m_i] << " \n";
                file << "# m_" << _lp.mass_labels[m_j] << "=" << format << _lp.masses[m_j]
                     << ", eps=" << format << _individual_naik_epsilons[m_j] << " \n";
                file << std::setw(0) << std::setprecision(0) << std::scientific;
                for (int ch = 1; ch <= 8; ch++){
                    for ( int momentum_index  = 0 ; momentum_index < _momentum_elements_per_axis[corr_axis_index]; momentum_index++ ){
                        for (int w = 0; w < _corr_ls[corr_axis_index]; w++){
                            LineFormatter output = file.tag(_lp.mass_labels[m_i]+_lp.mass_labels[m_j]+" M"+std::to_string(ch)
                                                            +" "+name_mapping[ch]);
                            output << w;
                            output << _momenta[ _momenta_start_index[corr_axis_index] + 3 * momentum_index ] ; 
                            output << _momenta[ _momenta_start_index[corr_axis_index] + 3 * momentum_index + 1 ] ;
                            output << _momenta[ _momenta_start_index[corr_axis_index] + 3 * momentum_index + 2 ] ;
                            output << real(_correlators[corr_index(w, ch, mass_index , corr_axis_index , momentum_index)]);
                            output << imag(_correlators[corr_index(w, ch, mass_index , corr_axis_index , momentum_index)]);
                        }
                    }
                }
                mass_index++;
            }
        }
        #undef format
    }
}

#define INITMEAS(floatT,LAYOUT,HALO,HALOSPIN,NSTACKS) \
template class measureHadrons<floatT, true, HALO, HALOSPIN, LAYOUT, NSTACKS, R18 , R14>;
INIT_PLHHSN(INITMEAS)

