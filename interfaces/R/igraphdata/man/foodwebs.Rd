\name{foodwebs}
\docType{data}
\alias{foodwebs}
\title{A collection of food webs}
\description{A list of graphs. Each one is a food web, i.e. a directed
  graph a predator-prey relationships.}
\usage{
foodwebs
}
\format{
  A named list of directed \code{igraph} graph objects. Here are the
  list of the graphs included:
  \describe{
    \item{\sQuote{ChesLower}}{Lower Chesapeake Bay in Summer.

      Reference: Hagy, J.D. (2002) Eutrophication, hypoxia and
      trophic transfer efficiency in Chesapeake Bay PhD Dissertation,
      University of Maryland at College Park (USA), 446 pp.}
    \item{\sQuote{ChesMiddle}}{Middle Chesapeake Bay in Summer.

      Reference: same as for \sQuote{ChesLower}.}
    \item{\sQuote{ChesUpper}}{Upper Chesapeake Bay in Summer.

      Reference: same as for \sQuote{ChesLower}.}
    \item{\sQuote{Chesapeake}}{Chesapeake Bay Mesohaline Network.

      Reference: Baird D. & Ulanowicz R.E. (1989) The seasonal dynamics
      of the Chesapeake Bay ecosystem. Ecological Monographs 59:329-364.}
    \item{\sQuote{CrystalC}}{Crystal River Creek (Control).

      Reference: Homer, M. and W.M. Kemp. Unpublished Ms. See also
      Ulanowicz, R.E. 1986. Growth and Development: Ecosystems
      Phenomenology. Springer, New York. pp 69-79.}
    \item{\sQuote{CrystalD}}{Crystal River Creek (Delta Temp).

      Reference: same as for \sQuote{CrystalD}.}
    \item{\sQuote{Maspalomas}}{Charca de Maspalomas.

      Reference: Almunia, J., G. Basterretxea, J. Aristegui, and
      R.E. Ulanowicz. (1999) Benthic- Pelagic switching in a coastal
      subtropical lagoon. Estuarine, Coastal and Shelf Science
      49:363-384.}
    \item{\sQuote{Michigan}}{Lake Michigan Control network.

      Reference: Krause, A. and D. Mason. (In preparation.) A. Krause,
      PhD. Dissertation, Michigan State University. Ann Arbor, MI. USA.}
    \item{\sQuote{Mondego}}{Mondego Estuary - Zostrea site.

      Reference: Patricio, J. (In Preparation) Master's
      Thesis. University of Coimbra, Coimbra, Portugal.}
    \item{\sQuote{Narragan}}{Narragansett Bay Model.

      Reference: Monaco, M.E. and R.E. Ulanowicz. (1997) Comparative
      ecosystem trophic structure of three U.S. Mid-Atlantic
      estuaries. Mar. Ecol. Prog. Ser. 161:239-254.}
    \item{\sQuote{Rhode}}{Rhode River Watershed - Water Budget.

      Reference: Correll, D. (Unpublished manuscript) Smithsonian
      Institute, Chesapeake Bay Center for Environmental Research,
      Edgewater, Maryland 21037-0028 USA.}
    \item{\sQuote{StMarks}}{St. Marks River (Florida) Flow network.

      Reference: Baird, D., J. Luczkovich and R. R. Christian. (1998)
      Assessment of spatial and temporal variability in ecosystem
      attributes of the St Marks National Wildlife Refuge, Apalachee Bay,
      Florida. Estuarine, Coastal, and Shelf Science 47: 329-349.}
    \item{\sQuote{baydry}}{Florida Bay Trophic Exchange Matrix, dry season.

      Reference: Ulanowicz, R. E., C. Bondavalli, and
      M. S. Egnotovich. 1998. Network analysis of trophic dynamics in
      South Florida ecosystems, FY 97: the Florida Bay ecosystem. Annual
      Report to the United States Geological Service Biological Resources
      Division, University of Miami Coral Gables, [UMCES] CBL 98-123,
      Maryland System Center for Environmental Science, Chesapeake
      Biological Laboratory, Maryland, USA.}
    \item{\sQuote{baywet}}{Florida Bay Trophic Exchange Matrix, wet season.
      
      Reference: same as for \sQuote{baydry}.}
    \item{\sQuote{cypdry}}{Cypress, dry season.

      Reference: Ulanowicz, R. E., C. Bondavalli, and
      M. S. Egnotovich. 1997. Network analysis of trophic dynamics in
      South Florida ecosystems, FY 96: the cypress wetland
      ecosystem. Annual Report to the United States Geological Service
      Biological Resources Division, University of Miami Coral Gables,
      [UM-CES] CBL 97-075, Maryland System Center for Environmental
      Science, Chesapeake Biological Laboratory.}
    \item{\sQuote{cypwet}}{Cypress, wet season.

      Reference: same as for \sQuote{cypdry}.}
    \item{\sQuote{gramdry}}{Everglades Graminoids - Dry Season.

      Reference: Ulanowicz, R. E., J. J. Heymans, and
      M. S. Egnotovich. 2000. Network analysis of trophic dynamics in
      South Florida ecosystems, FY 99: the graminoid ecosystem. Technical
      Report TS-191-99, Maryland System Center for Environmental Science,
      Chesapeake Biological Laboratory, Maryland, USA.}
    \item{\sQuote{gramwet}}{Everglades Graminoids - Wet Season.

      Reference: same as for \sQuote{gramdry}.}
    \item{\sQuote{mangdry}}{Mangrove Estuary, Dry Season.

      Reference: Ulanowicz, R. E., C. Bondavalli, J. J. Heymans, and
      M. S. Egnotovich. 1999. Network analysis of trophic dynamics in
      South Florida ecosystems, FY 98: the mangrove ecosystem. Technical
      Report TS-191-99, Maryland System Center for Environmental Science,
      Chesapeake Biological Laboratory, Maryland, USA.}
    \item{\sQuote{mangwet}}{Mangrove Estuary, Wet Season.
      
      Reference: same as for \sQuote{mangdry}.}

    Each graph has the following vertex attributes: \sQuote{name} is the
    name of the species, \sQuote{ECO} is the type of the node, and
    integer value between one and five, meaning:
    \enumerate{
      \item Living/producing compartment
      \item Other compartment
      \item Input
      \item Output
      \item Respiration.
    }
    The \sQuote{Biomass} vertex attribute contains the biomass of the
    species.

    Edges are weighted, and the weights denote energy flux between the
    species involved.

    The graphs also contain some informative graph attributes:
    \sQuote{Author}, \sQuote{Citation}, \sQuote{URL}, and
    \sQuote{name}.
  }
}
\source{See references for the individual webs above. The data itself
  was downloaded from
  \url{http://vlado.fmf.uni-lj.si/pub/networks/data/bio/foodweb/foodweb.htm}.
}
\references{See them above.}
\keyword{datasets}
