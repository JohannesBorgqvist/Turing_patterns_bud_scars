
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/function/Expression.h>
#include <dolfin/math/basic.h>
#include <Eigen/Dense>


// cmath functions
using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cosh;
using std::sinh;
using std::tanh;
using std::exp;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sqrt;
using std::ceil;
using std::fabs;
using std::floor;
using std::fmod;
using std::max;
using std::min;

const double pi = DOLFIN_PI;


namespace dolfin
{
  class dolfin_expression_cfb3d0893281eb42fb93ee5d628c6340 : public Expression
  {
     public:
       double ic_u;
double ic_v;
double sigma;


       dolfin_expression_cfb3d0893281eb42fb93ee5d628c6340()
       {
            _value_shape.push_back(2);
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          values[0] = ic_u+np.random.normal(0,sigma);
          values[1] = ic_v+np.random.normal(0,sigma);

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "ic_u") { ic_u = _value; return; }          if (name == "ic_v") { ic_v = _value; return; }          if (name == "sigma") { sigma = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "ic_u") return ic_u;          if (name == "ic_v") return ic_v;          if (name == "sigma") return sigma;
       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {

       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {

       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_cfb3d0893281eb42fb93ee5d628c6340()
{
  return new dolfin::dolfin_expression_cfb3d0893281eb42fb93ee5d628c6340;
}

