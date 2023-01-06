.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0( "\n", 
        paste( format( citation( "micEconCES" ), bibtex = FALSE ), collapse = "\n" ),
         "\nIf you have questions, suggestions, or comments ",
         "regarding the 'micEconCES' package, ",
         "please use the 'issue' tracker at the GitHub page of the package:\n",
         "https://github.com/micEcon/micEconCES" ),
      domain = NULL,  appendLF = TRUE )
}
