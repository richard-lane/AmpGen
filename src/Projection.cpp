#include "AmpGen/Projection.h"

#include "TH1D.h"
#include "TH2D.h"

using namespace AmpGen;

Projection::Projection( const std::function<double(const Event&)>& fcn, 
    const std::string& name, const std::string& xAxisTitle,
    const unsigned int& nBins, const double& min, const double& max,
    const std::string& units 
    ) :
  m_func(fcn),
  m_name(name),
  m_xAxisTitle(xAxisTitle),
  m_units(units),
  m_nBins(nBins), 
  m_min(min),
  m_max(max)
{
  m_width = (m_max-m_min)/double(m_nBins);
}

TH1D* Projection::plot(const std::string& prefix) const {
  std::string p = ( prefix==""?"":prefix+"_");
  TH1D* plot = new TH1D( (p + m_name ).c_str(),"",m_nBins,m_min,m_max);
  plot->GetXaxis()->SetTitle( m_xAxisTitle.c_str() );
  char buffer[100];
  if( m_units != "" ) sprintf(buffer,"\\mathrm{Entries} / (%0.2f %s)", m_width,m_units.c_str());
  else sprintf( buffer, "\\mathrm{Entries} / (%0.2f)",m_width);
  plot->GetYaxis()->SetTitle( buffer );
  plot->GetYaxis()->SetTitleOffset(1.35);
  plot->SetMarkerSize(0);
  plot->SetMinimum(0);
  return plot;
}
std::function<size_t(const Event& evt )> Projection::binFunctor() const {
  return [this](auto& evt){ 

    return int ( ( (*this)(evt) - m_min ) / m_width ) ; };
}

TH2D* Projection2D::plot(const std::string& prefix) const {
  std::string p = ( prefix==""?"":prefix+"_");
  TH2D* plot = new TH2D( ( p + xAxis.m_name +"_vs_"+yAxis.m_name).c_str(),"",
      xAxis.m_nBins,xAxis.m_min,xAxis.m_max ,
      yAxis.m_nBins,yAxis.m_min,yAxis.m_max );

  plot->GetXaxis()->SetTitle( xAxis.m_xAxisTitle.c_str() ); 
  plot->GetYaxis()->SetTitle( yAxis.m_xAxisTitle.c_str() ); 
  plot->GetYaxis()->SetTitleOffset(1.35);
  return plot; 
}

const std::string Projection::name() const { return m_name; }
double Projection::operator()( const Event& evt ) const { return m_func( evt ); }

std::pair<double, double> Projection2D::operator()( const Event& evt ) const
{
  return {xAxis.m_func( evt ), yAxis.m_func( evt )};
}
