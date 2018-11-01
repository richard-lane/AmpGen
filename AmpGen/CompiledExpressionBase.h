#ifndef AMPGEN_COMPILEDEXPRESSIONBASE_H
#define AMPGEN_COMPILEDEXPRESSIONBASE_H
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <future>

#include "AmpGen/Expression.h"
#include "AmpGen/CacheTransfer.h"

namespace AmpGen
{
  std::string programatic_name( std::string s );
  
  class MinuitParameter;
  class MinuitParameterSet;
  class ASTResolver; 

  /** @class CompiledExpressionBase
   *  Base class for compiled expressions, i.e. expressions that are (almost) ready to be evaluated.
   *  Handles (some) resolution and compilation behaviour, allows management of CompiledExpressions
   *  with explicitly referring to their return type, which is specified by template parameter
   *  in the implementation CompiledExpression
   */
  class CompiledExpressionBase
  {
  protected:
    Expression                                      m_obj;
    std::string                                     m_name;
    DebugSymbols                                    m_db;
    std::map<std::string, size_t>                   m_evtMap;
    std::shared_future<bool>*                       m_readyFlag;
    std::vector<std::pair<uint64_t, Expression>>    m_dependentSubexpressions;
    std::vector<std::shared_ptr<CacheTransfer>>     m_cacheTransfers;
  public:
    CompiledExpressionBase( const Expression& expression, 
                            const std::string& name,
                            const DebugSymbols& db=DebugSymbols(), 
                            const std::map<std::string,size_t>& evtMapping = {} );
    CompiledExpressionBase( const std::string& name );
    CompiledExpressionBase() : m_readyFlag(nullptr) {}

    void resolve(const MinuitParameterSet* mps = nullptr);
    void prepare();
    void compile(const std::string& fname="", const bool& wait = false ); 
    void to_stream( std::ostream& stream ) const;
    unsigned int hash() const;
    std::string name() const;
    virtual bool link( void* handle              )          = 0;
    virtual bool link( const std::string& handle )          = 0;
    virtual void setExternal( const double& value, const unsigned int& address ) = 0;
    virtual void resizeExternalCache( const size_t& N ) = 0;
    virtual bool isReady() const               = 0;
    virtual std::string returnTypename() const = 0;
    virtual std::string fcnSignature()   const = 0;
    virtual std::string args()           const = 0;
    virtual void print() const                 = 0;
    virtual ~CompiledExpressionBase() = default;
    virtual size_t returnTypeSize() const      = 0;    
  private:
    void addDebug( std::ostream& stream ) const;
    void addDependentExpressions( std::ostream& stream, size_t& sizeOfStream ) const;
    void resolveParameters( ASTResolver& resolver );
  };
  std::ostream& operator<<( std::ostream& os, const CompiledExpressionBase& expression );
}// namespace AmpGen
#endif