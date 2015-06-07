#include "include/Base.h"
#include "include/NtupleProducer.h"

ClassImp(Base)

Base::Base()
{
}

Base::~Base()
{
}

// first - parent idx, second - if after radiation
/*std::pair<int,bool> SKY::Base::truthGetParent(int p)
{
   if( p == -666 ) return std::make_pair(-666,false);
   
   int pPdgId = ntP->mc_pdgId->at(p);
      
   bool isDuplicate = true;
   bool afterRadiation = false;
   
   int parentIdx = -666;
   int parentPDGID = -666;
   
   while( isDuplicate )
     {	
	int parent = -666;
	int parentPdgId = -666;
	int childN = -666;
	
	if( ntP->mc_parent_index->at(p).size() > 0 )
	  {	
	     parent = ntP->mc_parent_index->at(p).at(0);
	     childN = ntP->mc_child_index->at(parent).size();
	     parentPdgId = ntP->mc_pdgId->at(parent);
	  }   
	
	isDuplicate = ( parentPdgId == pPdgId &&
			parent != -666 &&
			childN == 1 );

	if( childN == 2 && parent != -666 )
	  afterRadiation = ( parentPdgId == pPdgId &&
			     (ntP->mc_child_index->at(parent).at(0) == 22 ||
			      ntP->mc_child_index->at(parent).at(1) == 22) );
	
	if( afterRadiation && parent != -666 )
	  {
	     isDuplicate = true;
	  }
	
	parentIdx = parent;
	
	if( parentIdx >= 0 && ntP->mc_pdgId->size() > 0 )
	  parentPDGID = ntP->mc_pdgId->at(parentIdx);
	
	p = parent;
     }   
   
   return std::make_pair(parentIdx,afterRadiation);
}*/

float Base::GetDPhi(float phi1,float phi2)
{
   float deltaPhi=fabs(phi1-phi2);
   if (deltaPhi>M_PI) deltaPhi=2*M_PI-deltaPhi;
   return deltaPhi;
}

float Base::GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
   float DeltaPhi = TMath::Abs(phi2 - phi1);
      if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
