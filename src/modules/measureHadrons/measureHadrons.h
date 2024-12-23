//
// Created by luis on 11.01.21.
//

#include "../../gauge/gaugefield.h"
#include "../../base/latticeContainer.h"
#include "../hisq/hisqSmearing.h"
#include "../dslash/dslash.h"

#pragma once
enum source_type {
    POINT
};
enum action_type {
    HISQ , stdStag
};
enum correlator_axis{
    x_axis, y_axis, z_axis, t_axis
};

///convert string to enum for switch()
static std::map<std::string, source_type> source_type_map = {
        {std::string("POINT"), POINT}
};
static std::map<std::string, action_type> action_type_map = {
        {std::string("HISQ"), HISQ},
        {std::string("stdStag"), stdStag}
};
static std::map<std::string, correlator_axis> correlator_axis_map = {
        {std::string("x"), x_axis},
        {std::string("y"), y_axis},
        {std::string("z"), z_axis},
        {std::string("t"), t_axis}
};

template<class floatT>
struct measureHadronsParam : LatticeParameters {

    Parameter<std::string> action; //FIXME this feature isn't implemented yet-
    Parameter<std::string> source_type; //FIXME this feature isn't implemented yet

    DynamicParameter<std::string> correlator_axes;
    DynamicParameter<int> upper_momenta_per_correlator_axis;
    Parameter<std::string> measurements_dir;
    DynamicParameter<floatT> masses;
    DynamicParameter<std::string> mass_labels;

    Parameter<int, 4> source_coords; //x y z t
    Parameter<floatT> cg_residue;
    Parameter<int> cg_max_iter;

    //! you have to provide the same number of parameters as you did for the masses
    //! for each of the following (but they are optional)
    DynamicParameter<floatT> naik_epsilons_individual;
    DynamicParameter<floatT> cg_residues_individual;
    DynamicParameter<int> cg_max_iters_individual;

    Parameter<bool> use_unit_conf;

    // constructor
    measureHadronsParam() {
        add(measurements_dir, "measurements_dir");
        add(action, "action");
        add(masses, "masses");
        add(mass_labels, "mass_labels");
        add(source_type, "source_type");
        add(source_coords, "source_coords");
        add(correlator_axes, "correlator_axes");
        add(upper_momenta_per_correlator_axis, "upper_momenta_per_correlator_axis");

        addDefault(cg_residue, "cg_residue", static_cast<floatT>(1e-6));
        addDefault(cg_max_iter, "cg_max_iter", static_cast<int>(10000));

        addOptional(naik_epsilons_individual, "naik_epsilons_individual");
        addOptional(cg_residues_individual, "cg_residues_individual");
        addOptional(cg_max_iters_individual, "cg_max_iters_individual");

        addDefault(use_unit_conf, "use_unit_conf", false);
    }

    /*! call this after parameter read in */
    void check_for_nonsense() {
        if (masses.numberValues() != mass_labels.numberValues()){
            throw std::runtime_error(stdLogger.fatal("Number of mass labels is not equal to number of masses"));
        }
        if (naik_epsilons_individual.isSet() && (naik_epsilons_individual.numberValues() != masses.numberValues())){
            throw std::runtime_error(stdLogger.fatal("Number of naik_epsilons_individual is not equal to number of masses. "
                           "Use -1 as placeholder for no naik term."));
        }
        if (cg_residues_individual.isSet() && (cg_residues_individual.numberValues() != masses.numberValues())){
            throw std::runtime_error(stdLogger.fatal("Number of cg_residues_individual is not equal to number of masses. "
                           "Use -1 as placeholder for default value."));
        }
        if (cg_max_iters_individual.isSet() && (cg_max_iters_individual.numberValues() != masses.numberValues())){
            throw std::runtime_error(stdLogger.fatal("Number of cg_max_ites_individual is not equal to number of masses. "
                           "Use -1 as placeholder for default value."));
        }
        for (size_t i = 0; i < 4; i++){
            if (source_coords()[i] >= latDim()[i]){
                throw std::runtime_error(stdLogger.fatal("source_coords[", i, "] is greater than lattice extent!"));
            }
        }
        if ((source_coords()[0]+source_coords()[1]+source_coords()[2]+source_coords()[3]) % 2 != 0 ){
            throw std::runtime_error(stdLogger.fatal("Pointsource is odd but needs to be even!"));
        }
        for ( size_t i = 0 ; i <  correlator_axes.numberValues() ; i++ ){
            std::string correlator_axis = correlator_axes.get()[i] ;
            if (nodeDim()[correlator_axis_map[correlator_axis]] != 1){
                throw std::runtime_error(stdLogger.fatal("Don't split the lattice along the correlator axis!"));
            }
        }
    }
};


template<class floatT, bool onDevice, size_t HaloDepth, size_t HaloDepthSpin, Layout Source,  size_t NStacks, CompressionType compHISQ , CompressionType compstdStag>
class measureHadrons {

private:

    //Gaugefield<floatT, onDevice, HaloDepth, R18>& _gauge;
    CommunicationBase& _commBase;
    measureHadronsParam<floatT>& _lp;

    //! parameters derived from the ones set in the parameter struct:
    const size_t _n_masses;
    const size_t _n_channels = 8;
    const size_t _n_correlator_axes;
    const size_t _n_upper_momenta_per_axis;
    size_t _n_correlator_elements ;
    size_t _n_momenta_elements ;
    const source_type _src_type; //FIXME support for different sources
    const action_type _act_type; //FIXME support for wilson,clover fermions

    const bool _use_individual_naik_epsilons;
    const bool _use_individual_cg_max_iters;
    const bool _use_individual_cg_residues;
    std::vector<floatT> _individual_cg_residues;
    std::vector<int> _individual_cg_max_iters;
    std::vector<floatT> _individual_naik_epsilons;

    const size_t _vol4;
    std::array<size_t,4> _lat_extents; // "lx, ly, lz, lt"
    std::array<size_t,4> _global_lat_extents; // "LX, LY, LZ, LT"
    std::array<size_t,4> _tmp_axis_indices;
    std::vector<size_t> _axes_indices; // later of size 4*_n_correlator_axes where [4*_n_correlator_axes + 0] is corr axis, 4*_n_correlator_axes + (1 to 3)] are not corr axis
    std::vector<size_t> _corr_axes_start_index ;
    std::vector<size_t> _momenta_start_index ;
    std::vector<size_t> _momentum_elements_per_axis ;
    std::vector<size_t> _momenta ;

    std::vector<size_t> _corr_ls ;

    sitexyzt _current_source;

    //! In the end this contains the *contracted* propagators for each mass combination and spacetime point (vol4)
    std::vector<LatticeContainer<false,GPUcomplex<floatT>>> _contracted_propagators;

    //! In the end this contains the values of the correlator at each z for each mass combination and channel
    //! memory layout of this vector: see constructor and corr_index()
    std::vector<GPUcomplex<floatT>> _correlators;
    template <const size_t n_color> 
    void setup_memory_StaggeredType(std::vector<Spinorfield<floatT,false,All,0,n_color>> & prop_containers) ;
    template <class StaggeredTypeAction , const size_t n_color>
    void invert_StaggeredType_Dslash(StaggeredTypeAction & dslash , size_t massIndex , floatT mass , std::vector<Spinorfield<floatT,false,All,0,n_color>> & prop_containers) ;
    void compute_StaggeredType_correlator() ;
    template <const size_t n_color>
    void contract_prop_containers(std::vector<Spinorfield<floatT,false,All,0,n_color>> & prop_containers) ;

public:

    explicit measureHadrons(CommunicationBase& commBase, measureHadronsParam<floatT>& lp ) : //, Gaugefield<floatT,onDevice,HaloDepth,comp>& gauge) :
            _commBase(commBase),
            _lp(lp),
            _n_masses(_lp.masses.get().size()),
            _n_correlator_axes(_lp.correlator_axes.get().size()),
            _n_upper_momenta_per_axis(_lp.upper_momenta_per_correlator_axis.get().size()),
            _src_type(source_type_map[_lp.source_type()]),
            _act_type(action_type_map[_lp.action()]),
            _use_individual_naik_epsilons(lp.naik_epsilons_individual.isSet()),
            _use_individual_cg_max_iters(lp.cg_max_iters_individual.isSet()),
            _use_individual_cg_residues(lp.cg_residues_individual.isSet()),
            _individual_cg_residues(lp.cg_residues_individual.get()),
            _individual_cg_max_iters(lp.cg_max_iters_individual.get()),
            _individual_naik_epsilons(lp.naik_epsilons_individual.get()),
            _vol4(GIndexer<All,HaloDepth>::getLatData().vol4),
            _lat_extents({GIndexer<All,HaloDepth>::getLatData().lx,
                          GIndexer<All,HaloDepth>::getLatData().ly,
                          GIndexer<All,HaloDepth>::getLatData().lz,
                          GIndexer<All,HaloDepth>::getLatData().lt}),
            _global_lat_extents({GIndexer<All,HaloDepth>::getLatData().globLX,
                                GIndexer<All,HaloDepth>::getLatData().globLY,
                                GIndexer<All,HaloDepth>::getLatData().globLZ,
                                GIndexer<All,HaloDepth>::getLatData().globLT}),
            //_gauge(gauge),
            _current_source(lp.source_coords()[0], lp.source_coords()[1],lp.source_coords()[2],lp.source_coords()[3])
    {
        //! set up things which depend on _n_correlator_axes
        //! i.e. _axes_indices, _corr_ls, _n_correlator_elements
        _axes_indices.resize(4 * _n_correlator_axes);
        std::fill(_axes_indices.begin(), _axes_indices.end(), 0);
 
        _corr_ls.resize(_n_correlator_axes);
        _momentum_elements_per_axis.resize(_n_correlator_axes) ;
        _corr_axes_start_index.resize(_n_correlator_axes);
        _momenta_start_index.resize(_n_correlator_axes) ;
        _n_correlator_elements = 0 ;
        _n_momenta_elements    = 0 ;
        for (size_t i = 0; i < _n_correlator_axes; i++)
        {
            correlator_axis _tmp_corr_axis(correlator_axis_map[lp.correlator_axes.get()[i]]);
            switch(_tmp_corr_axis){
            case x_axis:
                _tmp_axis_indices = {0, 1, 2, 3};
                break;
            case y_axis:
                _tmp_axis_indices = {1, 0, 2, 3};
                break;
            case z_axis:
                _tmp_axis_indices = {2, 0, 1, 3};
                break;
            case t_axis:
                _tmp_axis_indices = {3, 0, 1, 2};
                break;
            default:
                throw std::runtime_error(stdLogger.fatal("Unknown correlator axis! Choose one from {x, y, z, t}."));
        }
            for (size_t j = 0; j < 4; j++)
            {
                _axes_indices[4 * i + j] = _tmp_axis_indices[j] ;
            }
            
            _corr_ls[i] = _lat_extents[_axes_indices[ i * 4 ]] ;
            _momentum_elements_per_axis[i] = ( _lp.upper_momenta_per_correlator_axis.get()[i] + 1 )
                                        * ( _lp.upper_momenta_per_correlator_axis.get()[i] + 1 )
                                        * ( _lp.upper_momenta_per_correlator_axis.get()[i] + 1 ) ;
            _n_momenta_elements    +=  3 * _momentum_elements_per_axis[i] ;
            _n_correlator_elements += _n_masses 
                                    * _n_masses 
                                    * _corr_ls[i] 
                                    * _n_channels
                                    * _momentum_elements_per_axis[i] ;
            
            if ( i == 0)
            {
                _corr_axes_start_index[i] = 0 ;
                _momenta_start_index[i]   = 0 ;
            }
            else
            {
                _corr_axes_start_index[i] = _n_masses * _n_masses * _corr_ls[i-1] * _n_channels
                                          * _momentum_elements_per_axis[ i - 1 ];
                _momenta_start_index[i]   = 3 * _momentum_elements_per_axis[ i - 1 ] ;
            }
            
            
        }
                
        

        //! set up _contracted_propagators
        for(size_t i = 0; i < _n_masses; i++) {
            for(size_t j = 0; j <= i; j++) {
                //! Here we need to give names without the "SHARED_" prefix to the MemoryManagement, otherwise all
                //! of these will point to the same memory.
                _contracted_propagators.emplace_back(std::move(LatticeContainer<false, GPUcomplex<floatT>>(_commBase,
                        "propmemA", "propmemB", "propmemC", "propmemD")));
                _contracted_propagators.back().adjustSize(_vol4);
            }
        }
        
        //! set up momenta vector
        _momenta.resize(_n_momenta_elements);
        std::fill(_momenta.begin(), _momenta.end(), 0);
        for (size_t i = 0; i < _n_upper_momenta_per_axis ; i++)
        {
            size_t momentum_index = 0 ;
            for (size_t k1 = 0; k1 <= _lp.upper_momenta_per_correlator_axis.get()[i] ; k1++)
            {
                for (size_t k2 = 0; k2 <= _lp.upper_momenta_per_correlator_axis.get()[i] ; k2++)
                {
                    for (size_t k3 = 0; k3 <= _lp.upper_momenta_per_correlator_axis.get()[i] ; k3++)
                    {
                    
                    
                        _momenta[_momenta_start_index[i] + 3 * momentum_index ]     = k1 ;
                        _momenta[_momenta_start_index[i] + 3 * momentum_index + 1 ] = k2 ;
                        _momenta[_momenta_start_index[i] + 3 * momentum_index + 2 ] = k3 ;

                        momentum_index++ ;
                    }
                }
            }
        }
        
        rootLogger.info( _n_momenta_elements , " " , _momenta.size() )  ;
        for (size_t i = 0; i < _n_momenta_elements ; i += 3)
        {
            rootLogger.info( _momenta[i] , _momenta[i + 1] , _momenta[i + 2] ) ;   
        }
        

        //! set up _correlators
        _correlators.resize(_n_correlator_elements);
        std::fill(_correlators.begin(), _correlators.end(), 0);

        //! set up source
        setCurrentSource(lp.source_coords());

        //! set up cg parameters
        if (_use_individual_cg_residues){
            for (floatT& val : _individual_cg_residues ){
                if (val < 0){ val = _lp.cg_residue(); } //! set to "default" value
            }
        } else {
            _individual_cg_residues = std::vector<floatT>(_n_masses,_lp.cg_residue());
        }

        if (_use_individual_cg_max_iters){
            for (int& val : _individual_cg_max_iters){
                if (val < 0){ val = _lp.cg_max_iter(); } //! set to "default" value
            }
        } else {
            _individual_cg_max_iters = std::vector<int>(_n_masses,_lp.cg_max_iter());
        }

        if (not _use_individual_naik_epsilons ){
            _individual_naik_epsilons = std::vector<floatT>(_n_masses,0.);
        }

    }

    void setCurrentSource(const int source[4]){
        _current_source.x = source[0];
        _current_source.y = source[1];
        _current_source.z = source[2];
        _current_source.t = source[3];
    }

    void compute_HISQ_correlators(Gaugefield<floatT,onDevice,HaloDepth,compHISQ>& gauge);
    void compute_stdStag_correlators(Gaugefield<floatT,onDevice,HaloDepth,compstdStag>& gauge);

    void write_correlators_to_file();

private:
    [[nodiscard]] inline int corr_index(const int w, const int channel, const int mass_index , const int corr_axis_index , const int momentum_index ) const{
        //! Channels start at 1 but arrays at 0, so it's channel-1
        return _corr_axes_start_index[corr_axis_index] 
                + mass_index * _n_channels * _momentum_elements_per_axis[corr_axis_index]  * _corr_ls[corr_axis_index] 
                + (channel - 1) * _momentum_elements_per_axis[corr_axis_index] * _corr_ls[corr_axis_index]   
                + momentum_index * _corr_ls[corr_axis_index]
                + w;
    }

    //FIXME implement arbitrary correlator axis support
    inline int phase_factor(const int channel, const int r, const int u, const int v){
        switch(channel){
            case 1:
                return 1 - 2*((r + u + v)%2);
            case 2:
                return 1;
            case 3:
                return 1 - 2*((u + v)%2);
            case 4:
                return 1 - 2*((r + v)%2);
            case 5:
                return 1 - 2*((r + u)%2);
            case 6:
                return 1 - 2*(r%2);
            case 7:
                return 1 - 2*(u%2);
            case 8:
                return 1 - 2*(v%2);
            default:
                throw std::runtime_error(stdLogger.fatal("Error in staggered phase factors: No such channel"));
        }
    }
};


