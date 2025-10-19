
#include "G4NCrystal/G4NCrystal.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManagerFactory.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysListFactory.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4PVPlacement.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4StepPoint.hh"
#include "G4ios.hh"
#include "G4GenericBiasingPhysics.hh"
#include "G4UserEventAction.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RandomDirection.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <optional>
#include <cassert>
#include <cstdlib> // for std::strtoull
#include "G4VisManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"  // Include the concrete implementation

#include "NCrystal/internal/utils/NCString.hh"//json output

//For convenience we simply reuse the lightweight hists from NCrystal (despite
//NCrystal issue #300 they are still handy):

#include "NCrystal/internal/utils/NCMath.hh"//missing include in NCHists (NCrystal issue #300)
#include "NCrystal/internal/utils/NCHists.hh"
using Hist1D = NCrystal::Hists::Hist1D<NCrystal::Hists::AllowWeights::YES,//YES due to NCrystal issue #300.
                                       NCrystal::Hists::OverflowHandling::Clamp//Clamp not Record record due to NCrystal issue #300
                                       >;

namespace {

  struct Arguments {
    bool do_help = false;
    bool do_vis = false;
    std::string shape = "box";
    double det_size_mm = 10.0;
    double film_thickness_um = 5.0;//5 allows us to cover the entire interesting
                                   //range (ends at ~4.2)
    double b10_enrichment = 0.95;
    unsigned long long nevents = 100000;
    static constexpr unsigned nbins = 400;
    unsigned long long seed = 123456;
    std::optional<std::string> outfile;
    bool isInteractive() const
    {
      return trailingArgs.has_value() || do_vis;
    }
    std::optional<std::vector<std::string>> trailingArgs;
    //derived convenience pars:
    std::string ncrystal_cfg_film;
    std::string ncrystal_cfg_gas;
  };

  void cfg2json( std::ostream& os, const Arguments& a )
  {
    namespace NC = NCrystal;
    NC::streamJSONDictEntry( os, "shape", a.shape, NC::JSONDictPos::FIRST );
    os << "\n    ";
    NC::streamJSONDictEntry( os, "det_size_mm", a.det_size_mm );
    os << "\n    ";
    NC::streamJSONDictEntry( os, "film_thickness_um", a.film_thickness_um );
    os << "\n    ";
    NC::streamJSONDictEntry( os, "b10_enrichment", a.b10_enrichment );
    os << "\n    ";
    NC::streamJSONDictEntry( os, "nevents", a.nevents );
    os << "\n    ";
    NC::streamJSONDictEntry( os, "nbins", a.nbins );
    os << "\n    ";
    NC::streamJSONDictEntry( os, "seed", a.seed );
    os << "\n    ";
    NC::streamJSONDictEntry( os, "ncrystal_cfg_film", a.ncrystal_cfg_film );
    os << "\n    ";
    NC::streamJSONDictEntry( os, "ncrystal_cfg_gas", a.ncrystal_cfg_gas,
                             NC::JSONDictPos::LAST );
  }

  std::vector<std::string> vis_mac() {
    return {
      "/tracking/verbose 2",
      "/vis/open",

      "/vis/viewer/set/background 0.1 0.1 0.1",

      "/vis/set/lineWidth 2",

      // # Disable auto refresh and quieten vis messages whilst scene and
      // # trajectories are established:
      "/vis/viewer/set/autoRefresh false",
      "/vis/verbose errors",
      // #
      // # Draw geometry:
      "/vis/drawVolume",
      // #
      // # Specify view angle:
      "/vis/viewer/set/viewpointVector -1 0 0",
      "/vis/viewer/set/lightsVector -1 0 0",
      // #
      // # Specify style (surface, wireframe, auxiliary edges,...)
      "/vis/viewer/set/style wireframe",
      "/vis/viewer/set/auxiliaryEdge true",
      "/vis/viewer/set/lineSegmentsPerCircle 60",
      // #
      // # Draw smooth trajectories at end of event, showing trajectory points
      // # as markers 2 pixels wide:
      "/vis/scene/add/trajectories smooth",
      "/vis/modeling/trajectories/create/drawByParticleID",
      "/vis/modeling/trajectories/drawByParticleID-0/set neutron green",
      "/vis/modeling/trajectories/drawByParticleID-0/set gamma yellow",
      "/vis/modeling/trajectories/drawByParticleID-0/set alpha blue",
      "/vis/modeling/trajectories/drawByParticleID-0/set Li7 red",
      "/vis/modeling/trajectories/drawByParticleID-0/set e- grey",
      "/vis/modeling/trajectories/drawByParticleID-0/set e+ grey",


      //      "/vis/modeling/trajectories/create/drawByCharge",
      //"/vis/modeling/trajectories/create/drawByParticleID",
      //"/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true",
      //"/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2",

      "/vis/scene/endOfEventAction accumulate",

      "/vis/scene/add/scale",         // # Simple scale line
      "/vis/scene/add/axes",          // # Simple axes: x=red, y=green, z=blue.

      "/vis/geometry/set/colour World -1 0 0 1",
      "/vis/geometry/set/colour film -1 1 0 0",
      "/vis/geometry/set/colour gas -1 0 1 0",
      //"/vis/viewer/set/style surface",
      //"/vis/viewer/set/hiddenMarker true",
      "/vis/viewer/set/viewpointThetaPhi 120 150",
      // #
      // # Re-establish auto refreshing and verbosity:
      "/vis/viewer/set/autoRefresh true",
      "/vis/verbose warnings",
    };
  }

  class GasSD final : public G4VSensitiveDetector {
  public:
    GasSD()
      : G4VSensitiveDetector("GasSD")
    {
    }

    void Initialize(G4HCofThisEvent* HCE) override
    {
      m_dep = 0.0;
      m_nsteps = 0;
    }

    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override
    {
      ++m_nsteps;
      auto dep = step->GetTotalEnergyDeposit();
      if (dep > 0)
        m_dep += dep;
      return true;
    }

    double getDep() const { return m_dep; }
    double getNSteps() const { return m_nsteps; }

  private:
    double m_dep = 0;
    unsigned long long m_nsteps = 0;
  };

  class FilmSD : public G4VSensitiveDetector {
    G4LogicalVolume * m_logVolGas = nullptr;
    void initLogVolPointer()
    {
      nc_assert_always(!m_logVolGas);
      auto lvs = G4LogicalVolumeStore::GetInstance();
      m_logVolGas = lvs->GetVolume("gas");
      nc_assert_always(m_logVolGas);
    }

  public:
    FilmSD()
      : G4VSensitiveDetector("FilmSD")
    {
    }

    void Initialize(G4HCofThisEvent* HCE) override
    {
      if (!m_logVolGas)
        initLogVolPointer();
      m_dist.reset();
      nc_assert_always( !m_dist.has_value() );
    }

    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override
    {
      if ( m_dist.has_value() )
        return false;
      auto track = step->GetTrack();
      if ( track->GetDynamicParticle()->GetPDGcode() != 2112 )
        return false;//only neutrons
      if ( track->GetTrackStatus() != fStopAndKill )
        return false;//is not being absorbed
      auto postStepPoint = step->GetPostStepPoint();
      nc_assert_always( postStepPoint->GetStepStatus() != fWorldBoundary );
      nc_assert_always( postStepPoint->GetStepStatus() != fGeomBoundary );

      G4ThreeVector p = postStepPoint->GetPosition();
      m_dist = m_logVolGas->GetSolid()->DistanceToIn(p);
      return true;
    }

    const std::optional<double>&
    getAbsnPtDist() const { return m_dist; }
  private:
    std::optional<double> m_dist;
  };

  class ThreadData {
    Hist1D hist_gasedep;//[keV]
    Hist1D hist_absptdist;//[micron]
    Hist1D hist_absptdist_detected;//[micron]
  public:
    ThreadData( const Arguments& args )
      : hist_gasedep( args.nbins, 0.0, 2000.0, "Gas Edep in keV"),
        hist_absptdist( args.nbins, 0.0, args.film_thickness_um * 1.2,
                        "Distance from absn pt to gas in micron"),
        hist_absptdist_detected( args.nbins, 0.0, args.film_thickness_um * 1.2,
                                "Distance from absn pt to gas in micron for"
                                " events exceeding Edep threshold")
    {
    }

    void addEvent( double edep, double absptdist )
    {
      hist_gasedep.fill(edep / CLHEP::keV);
      hist_absptdist.fill(absptdist / CLHEP::micrometer);
      if ( edep > 120*CLHEP::keV )
        hist_absptdist_detected.fill(absptdist / CLHEP::micrometer);
    }
    void merge(const ThreadData& o ) {
      hist_gasedep.merge(o.hist_gasedep);
      hist_absptdist.merge(o.hist_absptdist);
      hist_absptdist_detected.merge(o.hist_absptdist_detected);
    }

    void toJSON( std::ostream& os )
    {
      const bool include_data_arrays = true;

      os << "{ \"hists\" : { \n";

      os << "  \"gasedep\" : ";
      {
        std::ostringstream ss;
        hist_gasedep.toJSON( ss, include_data_arrays );
        os << ss.str();
        os << ",\n";
      }
      os << "  \"absptdist\" : ";
      {
        std::ostringstream ss;
        hist_absptdist.toJSON( ss, include_data_arrays );
        os << ss.str();
        os << "\n,";
      }
      os << "  \"absptdist_detected\" : ";
      {
        std::ostringstream ss;
        hist_absptdist_detected.toJSON( ss, include_data_arrays );
        os << ss.str();
        os << "\n";
      }
      os << "}\n }\n";
    }

  };

  struct ThreadDB {
    std::vector<std::shared_ptr<ThreadData>> threadDataList;
    std::mutex mtx;
  };
  ThreadDB& getThreadDB() {
    static ThreadDB db;
    return db;
  }
  std::shared_ptr<ThreadData> getNewThreadData( std::shared_ptr<Arguments> args ) {
    auto& db = getThreadDB();
    std::lock_guard<std::mutex> guard( db.mtx );
    auto nd = std::make_shared<ThreadData>( *args );
    db.threadDataList.emplace_back( nd );
    return nd;
  }

  std::shared_ptr<ThreadData> mergeAndFlushThreadData()
  {
    auto& db = getThreadDB();
    std::lock_guard<std::mutex> guard( db.mtx );
    auto it = db.threadDataList.begin();
    auto itE = db.threadDataList.end();
    if ( it == itE )
      throw std::runtime_error("mergeAndFlushThreadData: Nothing to merge");
    std::shared_ptr<ThreadData> res = *it;
    ++it;
    for ( ; it != itE; ++it )
      res->merge( **it );
    db.threadDataList.clear();
    return res;
  }
}

class MyEventAction final : public G4UserEventAction {
  std::shared_ptr<ThreadData> m_data;
  GasSD * m_gasSD = nullptr;
  FilmSD * m_filmSD = nullptr;
  void initSDPointers()
  {
    static std::mutex mtx;
    std::lock_guard<std::mutex> guard(mtx);
    nc_assert_always(!m_gasSD&&!m_filmSD);
    auto sdm = G4SDManager::GetSDMpointer();
    m_gasSD = dynamic_cast<GasSD*>(sdm->FindSensitiveDetector("GasSD"));
    m_filmSD = dynamic_cast<FilmSD*>(sdm->FindSensitiveDetector("FilmSD"));
    nc_assert_always(m_gasSD);
    nc_assert_always(m_filmSD);
  }
public:
  MyEventAction( std::shared_ptr<Arguments> args )
    : m_data( getNewThreadData( args  ) )
  {
  }

  void BeginOfEventAction(const G4Event*) override
  {
    if (!m_gasSD)
      initSDPointers();
    nc_assert_always(m_gasSD&&m_filmSD);


  }
  void EndOfEventAction(const G4Event* event) override {
    nc_assert_always(m_gasSD&&m_filmSD);

    auto absptdist = m_filmSD->getAbsnPtDist();
    if ( !absptdist.has_value() )
      return;//ignore events with no absorption

    const double edep = m_gasSD->getDep();
    m_data->addEvent( edep, absptdist.value() );
  }
};

class MyGeo final : public G4VUserDetectorConstruction {
  std::shared_ptr<Arguments> m_args;
public:

  MyGeo( std::shared_ptr<Arguments> args )
    : m_args( std::move(args) )
  {
  }

  G4VPhysicalVolume* Construct() override
  {
    G4Material * mat_vacuum = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Galactic", true );
    G4Material * mat_gas = G4NCrystal::createMaterial(m_args->ncrystal_cfg_gas);
    G4Material * mat_film = G4NCrystal::createMaterial(m_args->ncrystal_cfg_film);

    const double det_size = m_args->det_size_mm * CLHEP::mm;//half size
    const double film_thickness = m_args->film_thickness_um * CLHEP::micrometer;
    const double filmvol_size = det_size + film_thickness;
    const double d_world = 2.5 * (det_size + film_thickness);
    G4cout << "det_size[mm] = "<<det_size/CLHEP::mm<<G4endl;
    G4cout << "film_thickness[mm] = "<<film_thickness/CLHEP::mm<<G4endl;
    G4cout << "filmvol_size[mm] = "<<filmvol_size/CLHEP::mm<<G4endl;
    G4cout << "d_world[mm] = "<<d_world/CLHEP::mm<<G4endl;


    G4LogicalVolume* world_log = new G4LogicalVolume( new G4Box( "World", d_world, d_world, d_world ), mat_vacuum, "World", 0, 0, 0 );
    G4PVPlacement* world_phys = new G4PVPlacement( 0, G4ThreeVector(), world_log, "World", 0, false, 0 );
    G4LogicalVolume* film_log;
    G4LogicalVolume* gas_log;
    if ( m_args->shape == "sphere" ) {
      film_log = new G4LogicalVolume( new G4Sphere( "film", 0, filmvol_size,
                                                    0, CLHEP::twopi, 0, CLHEP::pi ), mat_film, "film", 0, 0, 0 );
      gas_log = new G4LogicalVolume( new G4Sphere( "gas", 0, det_size,
                                                   0, CLHEP::twopi, 0, CLHEP::pi ), mat_gas, "gas", 0, 0, 0 );
    } else if ( m_args->shape == "box" ) {
      film_log = new G4LogicalVolume( new G4Box( "film", filmvol_size, filmvol_size, filmvol_size ),
                                      mat_film, "film", 0, 0, 0 );
      gas_log = new G4LogicalVolume( new G4Box( "gas", det_size, det_size, det_size ),
                                     mat_gas, "gas", 0, 0, 0 );
    } else {
      throw std::runtime_error("unexpected shape");
    }
    new G4PVPlacement( 0, G4ThreeVector(), film_log, "film", world_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(), gas_log, "gas", film_log, false, 0 );
    m_logvols_ncrystal.emplace_back(film_log);
    m_logvols_ncrystal.emplace_back(gas_log);
    return world_phys;
  }
  void ConstructSDandField() override {
    auto sd_gas = new GasSD;
    auto sd_film = new FilmSD;
    G4SDManager::GetSDMpointer()->AddNewDetector( sd_gas );
    G4SDManager::GetSDMpointer()->AddNewDetector( sd_film );
    SetSensitiveDetector( "gas", sd_gas );
    SetSensitiveDetector( "film", sd_film );

    // Instantiate the biasing operator and attach to the volumes with NCrystal
    // materials (TODO: see if we can possibly automate this to simply bias any
    // volume where the material has an NCrystal scatter property:

    auto bias = new G4NCrystal::NCrystalBiasingOperator;
    for ( auto& e : m_logvols_ncrystal )
      bias->AttachTo( e );
  }
private:
  std::vector<G4LogicalVolume*> m_logvols_ncrystal;
};

class MyGun : public G4VUserPrimaryGeneratorAction {
  // Monochromatic source of neutrons, hitting the sample with initial direction (0, 0, 1)
  std::unique_ptr< G4ParticleGun > m_particleGun;
  double m_radius = 0.0;
public:
  MyGun( std::shared_ptr<Arguments> args )
    : m_particleGun( new G4ParticleGun(1) )
  {
    const double det_size = args->det_size_mm * CLHEP::mm;//half size
    const double film_thickness = args->film_thickness_um * CLHEP::micrometer;
    const double filmvol_size = det_size + film_thickness;
    const double d_world = 2.5 * (det_size + film_thickness);
    m_radius = 0.998*d_world;
    m_particleGun->SetParticleDefinition( G4ParticleTable::GetParticleTable()->FindParticle( "neutron" ) );
  }
  void GeneratePrimaries( G4Event* evt ) final
  {
    G4ThreeVector u = G4RandomDirection();
    G4ThreeVector pos;
    // if ( G4UniformRand()<0.5 )
    //   pos.set(0,0,0);
    // else
    pos = u*(-m_radius);
    const double neutron_wl = (G4UniformRand()*(10.0-0.5)+0.5)*CLHEP::angstrom;
    m_particleGun->SetParticleEnergy( neutronWavelengthToEKin( neutron_wl ) );
    m_particleGun->SetParticlePosition( pos );
    m_particleGun->SetParticleMomentumDirection( u );
    m_particleGun->GeneratePrimaryVertex( evt );
  }
  double neutronWavelengthToEKin( double wl ) {
    return 0.5 * CLHEP::h_Planck * CLHEP::h_Planck * CLHEP::c_squared / (wl*wl*CLHEP::neutron_mass_c2);
  }
private:
};

class MyActions final : public G4VUserActionInitialization {
  // Actions, for registering primary generator
  std::shared_ptr<Arguments> m_args;
public:
  MyActions( std::shared_ptr<Arguments> args )
    : m_args(args)
  {
  }
  void Build() const override
  {
    SetUserAction( new MyGun( m_args ) );
    SetUserAction( new MyEventAction( m_args ) );

  }
};

namespace {
  class CArgs final {
    int m_argc;
    char** m_argv;
  public:
    CArgs( const std::string& progname,
           const std::vector<std::string>& args )
      : m_argc(static_cast<int>(args.size()+1)),
        m_argv(new char*[m_argc + 1])//null terminated
    {
      for (int i = 0; i < m_argc+1; ++i)
        m_argv[i] = nullptr;
      m_argv[0] = new char[progname.size() + 1];
      std::strcpy(m_argv[0], progname.c_str());
      for (int i = 1; i < m_argc; ++i) {
        m_argv[i] = new char[args[i].size() + 1];
        std::strcpy(m_argv[i], args[i].c_str());
      }
    }

    ~CArgs() {
      for (int i = 0; i < m_argc; ++i)
        delete[] m_argv[i];
      delete[] m_argv;
    }

    int argc() const { return m_argc; }
    char** argv() { return m_argv; }//argv is not const

    //No copy/move:
    CArgs(const CArgs&) = delete;
    CArgs& operator=(const CArgs&) = delete;
    CArgs(const CArgs&&) = delete;
    CArgs& operator=(const CArgs&&) = delete;
  };


  void print_usage( std::ostream& os, const std::string& progname ) {
    os << "usage:\n\n";
    os << progname<< " [--h|--help] [--seed VAL] [-n VAL]";
    os << " [--out FILENAME]";
    os << " [--detsize-mm VAL]";
    os << " [--filmthickness-um VAL]";
    os << " [--b10enrich VAL]";
    os << " [--vis] [-- <arguments for G4UIExecutive>]";
    os << "\n";
  }

  Arguments parseArguments(int argc, char** argv) {
    Arguments args;
    for (int i = 1; i < argc; ++i) {
      const std::string a(argv[i]);
      if ( args.trailingArgs.has_value() ) {
        //Just consume trailing
        args.trailingArgs.value().push_back(argv[i]);
        continue;
      }
      if ( a == "--help" || a == "-h" ) {
        args.do_help = true;
      } else if ( a == "--vis" ) {
        //NB: use export QT_SCALE_FACTOR=2 for readable fonts on 4K displays
        args.do_vis = true;
      } else if (a == "--seed") {
        if (i + 1 >= argc)
          throw std::invalid_argument("--seed option requires a value.");
        args.seed = std::strtoull(argv[++i], nullptr, 10);
      } else if (a == "-n") {
        if (i + 1 >= argc)
          throw std::invalid_argument("-n option requires a value.");
        args.nevents = std::strtoull(argv[++i], nullptr, 10);
      } else if (a == "--out") {
        if (i + 1 >= argc)
          throw std::invalid_argument("--output option requires a value.");
        args.outfile = argv[++i];
      } else if (std::string(argv[i]) == "--") {
        //simply collect rest of args into trailingArgs list:
        args.trailingArgs.emplace();

      } else if (a == "--detsize-mm") {
        double val(-1.0);
        if (i + 1 < argc)
          val = std::stod(argv[++i]);
        if ( !(val > 0.0))
          throw std::invalid_argument("invalid or missing value for --detsize-mm option.");
        args.det_size_mm = val;
      } else if (a == "--filmthickness-um") {
        double val(-1.0);
        if (i + 1 < argc)
          val = std::stod(argv[++i]);
        if ( !(val > 0.0))
          throw std::invalid_argument("invalid or missing value for --filmthickness-um option.");
        args.film_thickness_um = val;
      } else if (a == "--b10enrich") {
        double val(-1.0);
        if (i + 1 < argc)
          val = std::stod(argv[++i]);
        if ( !(val > 0.0) || !(val<1.0))
          throw std::invalid_argument("invalid or missing value for --b10enrich option.");
        args.b10_enrichment = val;
      } else if (a == "--shape") {
        std::string val;
        if (i + 1 < argc)
          val = argv[++i];
        if ( !( val == "box" || val == "sphere" ) )
          throw std::invalid_argument("invalid or missing value for --shape option.");
        args.shape = val;
      } else {
        std::ostringstream ss;
        ss<<"Unsupported arguments: \""<<a<<"\"";
        throw std::invalid_argument(ss.str());
      }
    }

    std::string ncrystal_cfg_film;
    {
      std::ostringstream ss;
      ss << "stdlib::B4C_sg166_BoronCarbide.ncmat;temp=293.15K;atomdb=B is "
         << args.b10_enrichment<<" B10 "<<(1.0-args.b10_enrichment)<<" B11";
      args.ncrystal_cfg_film = ss.str();
    }
    args.ncrystal_cfg_gas = "gasmix::0.7xAr+0.3xCO2/massfractions/1.0atm/293.15K";

    return args;
  }
}

int main( int argc, char** argv ) {
  const auto args = std::make_shared<Arguments>(parseArguments(argc,argv));
  if ( args->do_help ) {
    print_usage(std::cout,argv[0]);
    return 0;
  }


  G4VisManager* vismgr = nullptr;
  G4UIExecutive* ui = nullptr;

  // We simply use the runmanager factory and use the env var
  // G4FORCE_RUN_MANAGER_TYPE if we want to switch MT/Serial/Tasking mode (and
  // G4FORCENUMBEROFTHREADS if needed):

  auto* runManager = G4RunManagerFactory::CreateRunManager( G4RunManagerType::Default );

  CLHEP::HepRandom::setTheSeed( args->seed );
  runManager->SetUserInitialization( new MyGeo(args) );  // Setup geometry

  // Setup HP physics-list (for NCrystal) with EMZ for more precise ionisation physics:
  G4VModularPhysicsList* physicsList = G4PhysListFactory().GetReferencePhysList( "QGSP_BIC_HP_EMZ" );

  //Ensure biasing physics is enabled for neutrons:
  G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics;
  biasingPhysics->Bias( "neutron" );
  physicsList->RegisterPhysics( biasingPhysics );
  runManager->SetUserInitialization( physicsList );

  runManager->SetUserInitialization( new MyActions( args ) );

  std::optional<CArgs> cargs;
  if ( args->isInteractive() ) {
    cargs.emplace( "myg4app",
                   args->trailingArgs.value_or(std::vector<std::string>{}) );
    ui = new G4UIExecutive( cargs.value().argc(), cargs.value().argv() );
    if (args->do_vis) {
      vismgr = new G4VisExecutive;
      vismgr->Initialize();
    }
  }

  runManager->Initialize();  // Initialize g4 run manager

  if ( args->isInteractive() ) {
    if ( args->do_vis  )
      for ( auto& e : vis_mac() )
        G4UImanager::GetUIpointer()->ApplyCommand(e);
    ui->SessionStart();
  } else {
    runManager->BeamOn( args->nevents );  // Perform simulations
  }
  delete vismgr;
  //FIXME??? delete ui;
  delete runManager;

  auto data = mergeAndFlushThreadData();

  std::string json;
  {
    std::ostringstream ss;
    ss << "{ \"simoutput\" : \n";
    data->toJSON( ss );
    ss << ",  \n";
    ss << "  \"cfg\" : \n";
    cfg2json( ss, *args );
    ss << " }\n";
    json = ss.str();
  }
  if ( args->outfile.has_value() ) {
    {
      std::ofstream outFile(args->outfile.value());
      if (!outFile)
        throw std::runtime_error("Could not open output file for writing");
      outFile << json;
    }
    std::cout<<"Wrote: "<<args->outfile.value()<<std::endl;
  } else {
    std::cout<<json;
  }

  return 0;
}
