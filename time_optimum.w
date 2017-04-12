\pdfoutput=1
\documentclass[11pt]{cweb}%
\usepackage[english]{babel}
\usepackage[T1]{fontenc}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{tikzsymbols}
\usepackage{fullpage}
\usepackage{bbold}%
\usepackage{bm}
%%\usepackage{CountriesOfEurope}
\usepackage{natbib}
\usepackage[all]{xy}%
\usepackage{color}%
\usepackage{rotating}%
\usepackage{a4wide,fullpage}%
\usepackage{setspace}%
\usepackage{enumerate}%
\usepackage{wasysym}
\usepackage{textcomp}
\usepackage{hyperref}
\setstretch{1.2}%
\newcommand{\one}[1]{\ensuremath{\mathbb{1}_{\left( #1 \right)}  } }%%
\newcommand{\g}{\,\boldsymbol{|}\,}%
\newcommand{\EE}[1]{\mathbb{E}\left[ #1 \right]}%
\newcommand{\pr}[1]{\ensuremath{\mathbb{P}\left( #1 \right) } }%%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\T}{\ensuremath{\mathcal{T}}}%
\newcommand{\be}{\begin{equation}}%
\newcommand{\ee}{\end{equation}}%
\newcommand{\tG}{\scalebox{1.2}{\tt G}}%
\newcommand{\tB}{\scalebox{1.2}{\tt B}}%
\newcommand{\tK}{\scalebox{1.2}{\tt K}}%
\newcommand{\norm}[2]{\ensuremath{\boldsymbol{|}{#1}\boldsymbol{|}_{#2} } }%
\newcommand{\In}{\ensuremath{\mathcal{I}_n} }%
\newcommand{\IN}{\ensuremath{\mathbb{N}}}%
\newcommand{\R}{\ensuremath{\mathbb{R}} }%
\newcommand{\bd}{\begin{displaymath}}%
\newcommand{\ed}{\end{displaymath}}%
\newcommand{\uZ}{\ensuremath{\underline{Z}}}%
\newcommand{\uR}{\ensuremath{\underline{R}}}%
\newcommand{\uX}{\ensuremath{\underline{X}}}%
\newcommand{\uY}{\ensuremath{\underline{Y}}}%
\newcommand{\cleanup}{ // clear all used memory: }%
\newcommand{\bone}[1]{\ensuremath{ \mathbb{1}\left( #1 \right) } }%
\newcommand{\EP}[2]{\ensuremath{\mathbb{E}^{#1}\left[ #2 \right] } }%
\newcommand{\hj}{\ensuremath{\hat{\jm}}}%
\title{\bf Viability selection and high fecundity\\ time to reach high frequency }
\author{ CWEB technical report\\ bjarki eldon\\ Museum f\"ur Naturkunde \\ 
{Leibniz Institut f\"ur Evolutions- und
Biodiversit\"atsforschung} \\  Berlin, Germany}
\date{\today }%
\begin{document}
\maketitle

\begin{abstract}
This code simulates viability selection  in a haploid population characterised
 by high fecundity and sweepstakes reproduction (HFSR).  An estimate of the
 expected  time of an allelic type to reach a given frequency, conditional on
 the event  of reaching high frequency. We exclude mutation. 
 This CWEB \citep{knuth1994cweb}
technical report describes corresponding C \citep{kernighan1988c} code.  CWEB documents may be compiled with {\tt cweave} and {\tt
ctangle}.
\end{abstract}


\tableofcontents



@q generate shasum: shasum -a 512224 file.tgz > file.tgz.sha1 @>
@q check shasum: shasum -a 512224 -c file.tgz.sha1 @>


@* {\bf Copyright}. 

Copyright {\copyright} \the\year{} Bjarki Eldon \newline

This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.




@* {\bf Introduction}. 

Some populations are highly fecund broadcast spawners  and may be characterised by Type III
survivorship curve.    The reproduction mode  of such populations has been
described as sweepstakes reproduction  where  few parents  contribute most of
the offspring  to a new generation.  Reproduction models which take into
account sweepstakes reproduction  do so through heavy-tailed, or skewed,
offspring distributions.   The impact of such  reproduction modes  on
selection has been little discussed \citep{Der2012,foucart2013impact,EGT10}.   


We consider a new model of  HFSR  in a haploid population of fixed size $N$.
 In each generation, 
 individual $i$ for $i\in [N] :=  \{1, 2, \ldots, N\}$ for   $N \in \IN := \{1,2,\ldots \}$   
independently  contributes  a  random number $X_i$ of  juveniles.   If the
total count of juveniles exceeds $N$  random sampling  of juveniles takes
place in which $N$ juveniles are sampled to form the new set of adults.      In case of a
highly fecund population  with sweepstakes reproduction (HFSR population), the distribution of
$X_i$  is heavy-tailed with
parameters  $\alpha,C, \gamma > 0$  and  mass function  
\be\label{eq:skew}
   \mathbb{P}\left( X =k \right) :=  C\left(  \frac{1}{k^\alpha}  -   \frac{
    1 }{ (k + 1)^\alpha } \right), \quad 1 \le k  \le     \gamma.
\ee
One can choose $C$ so that $\mathbb{P}(X_1 = 0) \ge 0$ and $\mathbb{E}[X_1] >
 1$.     Our main requirement is that $\mathbb{E}[X_1] > 1$  since then  the
 total number of juveniles is  at least $N$ with high probability for large
 $N$. 


We model viability selection as follows. We assume there are $n$ allelic types
segregating in the population; we label these types by the typespace $E =
\{0,1, \ldots, n-1\}$.  The juveniles inherit the types of their  parents
since we exclude mutation.     We assume there is a  \emph{trait
function} which maps the genetic type to a trait value.  We assume  the trait
function 
\be\label{eq:z} 
 z(i) = \frac{i}{i + 1}, \quad i \in \{0, 1, \ldots, n-1\}.
\ee
We assume there is a \emph{fitness function} which maps the  trait value to a
fitness value.  We consider an exponential fitness function, where $s$ denotes
the strength of selection and $z_0$ the optimal trait value, 
\be\label{eq:exp}
    w(z) =  \exp\left( - s(z - z_0)^2 \right), \quad z \in [0,1];
\ee
and  an algebraic fitness function
\be\label{eq:alg}
    w(z) =  \frac{1}{ 1 +  s(z - z_0) }, \quad z \in [0,1].
\ee
  If the count of juveniles   is greater than $N$  we draw  a random
  exponential with rate $w(z)$  from either the algebraic \eqref{eq:alg} or
  exponential \eqref{eq:exp} fitness function.   The $N$ juveniles with
  smallest  times  then form the new set of adults.  If the count of juveniles
  equals $N$ then all juveniles survive; we draw a new set of juveniles in
  case the count is less than $N$.   


  Let $Y_r$ denote the frequency of  the type conferring highest fitness at
  time (generation) $r$. Define $T := \inf \{ r \in \IN : Y_r \ge y \}$
  as the first time $Y_r$ is at least  $y$; ie.\ the fittest type has reached
  frequency $y$.   We   estimate the conditional expected  time $\tau := \mathbb{E}\left[ T \, :  \,
  Y_r > 0 \,\, \forall \,\, r  \right] $.      For comparison with our HFSR
  model \eqref{eq:skew}  we  model the number of juveniles according to a
  Poisson distribution with mean
  $\mathbb{E}^{(\textrm{HFSR})}\left[X_1\right]$.   

@* {\bf Compile and run }.  

Use {\tt cweave} on the {\tt .w} file to generate {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file.

The necessary parameters have preset values (see section~\ref{sec:main}).  By
way of example,  assuming the executable is {\tt a.out}, then with random seed
{\tt 12345} the  command 
\begin{verbatim}
./a.out 12345 out.out 
\end{verbatim}
writes into {\tt out.out}
\begin{verbatim}
24
10
22
15
15
31
6
9
14
10
\end{verbatim}
ten realisations of the time $T$.  
@* {\bf Code}. 



@*1 {\bf Random number generator }. 

A  random number generator of choice is declaired  using the
$GSL\_RNG\_TYPE$ environment variable.  The default generator is the
`Mersenne~Twister' random number generator \cite{MN1998} as
implemented in  GSL.


@<random number generator@>=@#

@t declare the random number generator $rngtype$ @>@#
gsl_rng * rngtype ; @#
@t Define the function $setup\_rng$ which initializes $rngtype$:@>@#
void setup_rng( unsigned long int seed ) @#
{@#
   @q *gsl_rng_alloc( const gsl_rng_type *T ) @>
  @q const gsl_rng_type *T ; @>
   
   @q T = gsl_rng_default ; @>
  @q rngtype = gsl_rng_alloc( T ) ; @>
  
  @t set the type as $mt19937$ @>@#
  rngtype = gsl_rng_alloc(gsl_rng_mt19937); @#
  gsl_rng_set(rngtype,  seed); 
  
      gsl_rng_env_setup();  @#
   @q gsl_rng_default_seed = rngseed ; @>
  @q  printf("seed %lu\n", gsl_rng_default_seed); @>
}



@*1 {\bf Definitions}. 

@<object definitions@>=@#
#define MAX_JUVENILES 10000000



@*1 {\bf Draw values for $X_i$ }. 

 Draw values for $X_i$; the diploid juveniles. We take $C, \alpha > 0$ and consider 
\bd
  \pr{X_i = k} =  Ck^{-\alpha}   -  C\bone{k < \psi}(1+k)^{-\alpha} , \quad k \in [\psi], 
\ed 
and we observe that $\pr{X_i = 0} = 1 - C$.  To have the mean  $\mathbb{E}[X_1] > 1$ we require approximately  
$C > \alpha - 1$.      

@<initialize distribution for $X_i$@>=@#
void drawXi( int N,   int psi, double a, double b, gsl_ran_discrete_t * Pmass,    int * tXi,   gsl_rng * r)
{@#

        int k, teljari ; @#
    tXi[0] = 0; @#
    teljari = 0 ; @#
  while( (tXi[0] < N)  &&  (teljari < 1000000) ){@#
  teljari = teljari + 1 ; @#
   tXi[0] = 0; @#
   for( k = 1 ; k <= N ; k++){@#
     tXi[k] = (a > 0. ? (int)gsl_ran_discrete( r, Pmass) : gsl_ran_poisson(r, b) ) ; @#
     tXi[0] = tXi[0]  +  tXi[k]; }} @#

 assert( teljari < 1000000); @#

     assert( tXi[0] <= MAX_JUVENILES); @#
}




@*1 {\bf Update population}. 

Update population given numbers   $x_i$ of juveniles generated by each
individual. 

@<population update@>=@#

double update_population( int N, double variance,  int nalleles,   double s,
double znull,  double epsilon, int * Pop,  int * tXi, int * tempJuve,   double
*Z, double *locuseffects,  double *etimes, size_t * aindex,   gsl_rng * r )
{@#

          @t \newline @>@#
          /*  $N$ is number of pairs,  $L$ is  number of loci,  $s$ is
          selection coefficient,  $znull$ is  trait  optimum */
            @t \newline @>@#


  int i,  k, xindex; @#
  double Zbar = 0. ; @#
  double w;@#

  @t \newline @>@#
  /*  $tXi[0] =  X_1 + \cdots +  X_N$ is the total number of juveniles */
    @t \newline @>@#
    xindex = 0 ; @#
  for( i = 1 ; i <= N ; i++){ @#
    @t \newline @>@#
    /*  check if  individual  $i$ produced potential offspring, ie.\ if $X_i >
  0$ */
      @t \newline @>@#
      
   if( tXi[i] > 0 ){ @#
    for( k = 1 ; k <= tXi[i] ; k++){ @#
          @t \newline @>@#
          /* if no mutation, copy the type of the parent  */
            @t \newline @>@#
          tempJuve[xindex + 1]  =  Pop[i] ; @#
            @t \newline @>@#
          /*  compute the trait value $z_{i} =  \bone{v > 0} G(0,v)  +
    \tfrac{i}{1+i}$ of juvenile  $i$ where $v$ denotes the variance    */
            @t \newline @>@#
          Z[xindex + 1]  =  ( variance > 0. ? gsl_ran_gaussian_ziggurat( r,
    variance) : 0.)   +   locuseffects[ tempJuve[xindex + 1] ]; @#
           @t \newline @>@#
           /*  compute fitness value $w = 1/(1 + s(z_i - z_0)^2) $;  now exponential  fitness function $w = \exp\left(- s(z_i - z_0)^2 \right)$ */
             @t \newline @>@#
            @q  w =  gsl_sf_exp( -s*  gsl_pow_2( Z[xindex+1]  -  znull) ); @>
              w =  1./( 1. +  (s * gsl_pow_2( Z[xindex+1]  -  znull)) ) ; @#
             assert( w > 0); @#

             /* draw exponential times with rate $w$ */
                @t \newline @>@#
             etimes[xindex]  =  gsl_ran_exponential(r,  1./w); @#
           xindex = xindex + 1 ; } } }
         assert ( xindex == tXi[0] ); @#
      @t \newline @>
      /* sort the exponential times, if $tXi[0] > N$  */
          @t \newline @>
      if ( tXi[0] > N){ @#
          gsl_sort_index( aindex,  etimes,  1,  tXi[0]); @#
          @t \newline @>
       /* the first $N$ indexes in $aindex$  are the indexes of the surviving juveniles */
          @t \newline @>
           for( i = 0 ; i < N ; i++){ @#
             Pop[i+1]  =  tempJuve[ aindex[i] ];  @#
               @t \newline @>
                /* the trait value of the population is given by $\overline{z}
    =  \tfrac 1N \sum_i z_{\sigma(i)}$ where $\sigma(i)$ is the ordered index
    $i$, in ascending order of the associated exponential times; we compute and return the fraction of the null type, the most fit type */
                 @t \newline @> 
             
             @q Zbar = Zbar  +  ( Z[ aindex[i] ] / ( (double)N ) );  @>
              @t \newline @> 
              /* $Zbar$ counts the number of  alleles of the fittest type; either type $n-1$ or $0$ can be the fittest types  */
               @t \newline @> 
             Zbar = Zbar  +  (znull > 0. ? (Pop[i+1] < nalleles - 1 ? 0.0 : 1.0) : (Pop[i+1] > 0 ? 0. : 1.0))/( (double)N ); @#
             }}
         else{  @#
           @t \newline @>
            /*  exactly $N$ juveniles, so all survive */
                @t \newline @>
                  for( i = 0 ; i < N ; i++){ @#
             Pop[i+1]  =  tempJuve[i+1];  @#
             
           @q Zbar = Zbar +  (Z[i+1]/( (double)N ));  @>
             @t \newline @> 
              /* $Zbar$ counts the number of  alleles of the fittest type; either type $n-1$ or $0$ can be the fittest types  */
               @t \newline @> 
               Zbar = Zbar  +   (znull > 0. ? (Pop[i+1] < nalleles - 1 ? 0.0 : 1.0) : (Pop[i+1] > 0 ? 0. : 1.0) )/( (double)N ) ; @#
                  }}

  return( Zbar ) ; @#
}




@*1 {\bf Simulator}. 


Run many replicates.  

@<replicates@>=@#
void   simulator( int N, int nalleles,  double variance,     double a, double b,  int Psi, double
s, double znull, double epsilon, int nruns,   char skra[200],   gsl_rng * r)
{ @#

@q   printf("pop size %d\n", N ); @>
@q  printf("alleles %d\n", nalleles ); @>
@q  printf("variance %g\n", variance ); @>
@q printf("a %g\n", a ); @>
@q printf("b %g\n", b ); @>
@q  printf(" psi %d\n", Psi ); @>
@q printf("s %g\n", s ); @>
@q printf("zn %g\n", znull ); @>
@q printf("eps %g\n", epsilon ); @>
@q printf("runs %d\n", nruns ); @>


 @t \newline @>
    /* $N$ is population size; $nalleles$ is number of alleles */
       @t \newline @>

       double zbar ; @#

   int * Pop  =  (int * )calloc( N + 1, sizeof( int)); @#
    double * Z =   (double *)calloc(MAX_JUVENILES, sizeof(double)); @#
    
   size_t * aindex  =  (size_t *)calloc(MAX_JUVENILES, sizeof(size_t)); @#
   double * etimes =   (double *)calloc(MAX_JUVENILES, sizeof(double)); @#
  int * tempJuve = (int *)calloc( MAX_JUVENILES, sizeof(int)); @#
   double * PXi =  (double *)calloc( 1 + Psi, sizeof(double)); @#
  double * leffects  =  (double *)calloc( nalleles,  sizeof(double)); @#
  double X0; @#

    int * tXi  =  (int * )calloc( N + 1, sizeof( int)); @#
       
   int k, ngens ; @#
  double mean = 0. ; @#
  PXi[0] = 0. ; @#
  for( k = 1 ; k <=  Psi ; k++){ @#
	@t \newline @>
	/* $P(X_i = k) =  k^{-\alpha} - (k+1)^{-\alpha}$ for $1 \le k \le \psi$ */
	@t \newline @>
    PXi[k] =  (a > 0. ? (  pow( 1./( ((double)k) ), a) -  pow( 1./( ((double)(1+k))  ), a)) : 1. ) ; @#

    assert( PXi[k] >= 0.) ; @#

    mean = mean  +  ( ((double)k) * PXi[k] ); } @#
    
   @q   assert( mean >= 1. ); @>

   gsl_ran_discrete_t *  Pmass = gsl_ran_discrete_preproc( 1 + Psi, PXi); @#
   
   @t \newline @>
   /* here we set the allelic  type   effects $\xi_j$; one option might be $\xi_j = j/(1
   + j)$, another option might be $\xi_j =  j/n$ where $n$ is the number of
   types, and $0 \le j \le n-1$.  */
    @t \newline @>
   for( k = 0 ; k < nalleles ; k++){ @#
     leffects[k]  =  ((double)k) / ( (double)( 1 + k)) ; } @#

  int rep = 0; @#

  while( rep < nruns){ @#
  zbar = 0. ; @#
   @t \newline @>
    /* initialise population by assigning allelic type from $\{0, 1, \ldots,
  n-1\}$, where $n$ is number of types,   uniformly at random to
  each individual.  Initialise    $zbar  =  \overline{z} =   \tfrac 1N \sum_i
  z_i$ where $z_i =  \bone{\sigma >0}N(0,\sigma) +  \xi_i( g_i)$ where $g_i$
  is the genotype of individual $i$, and $N(0, \sigma)$ is a random Gaussian
  with mean 0 and variance $\sigma$  */
       @t \newline @>
       X0 = 0.0 ; @#
       for( k = 1 ; k <= N ; k++){ @#
        @t \newline @>
        /* assign a type modulo $n$ where $n$ is  number of types; set $\mathbb{R}_n := \{0, 1, \ldots, n-1\}$ and we assign type $a_j =  j \mod n$ where $a_j \in \mathbb{R}_n  $  */
  @t \newline @>
           Pop[k]  =  (int)(k%nalleles) ; @#
        @t \newline @>
             /* starts almost fixed at  type $n-1$; one copy of each of other alleles  */
            @t \newline @>
            @q  Pop[k]  =  (k-1 < nalleles ? k-1 : nalleles - 1); @>

           assert( Pop[k]  >=  0 ) ; @#
           assert( Pop[k] < nalleles) ; @#
          
          X0 = X0  +   (znull > 0. ? (Pop[k] < nalleles - 1 ? 0.0 : 1.0) : (Pop[k] > 0 ? 0. : 1.))/( (double)N ) ; @#
 
        zbar = zbar  +   ( variance > 0. ? gsl_ran_gaussian_ziggurat( r,
        variance) : 0.)  +  leffects[ Pop[k] ] ; @#

           } @#

         ngens  =  0; @#
             @t \newline @> 

         @t \newline @>
         /* $\varepsilon$ is  the fraction of  the null type - the most fit type */
           @t \newline @>

         while( ((ngens < 100000) && ( X0 < epsilon )) && (X0 > 0.0 ) ){ @#
          
         @q  printf("%g\n", zbar) ; @>

          drawXi(N, Psi, a, b, Pmass, tXi, r); @#

          zbar  =  update_population( N, variance,  nalleles,   s,  znull,  epsilon, Pop,  tXi,
          tempJuve, Z, leffects, etimes, aindex, r); @#
          ngens = ngens  +  1; @#
          X0 = zbar ; @#
          } @#
          
          if( X0 > 0.0){ @#
          FILE * f = fopen(skra, "a"); @#
           fprintf(f,  "%d\n", (X0 > 0 ? ngens : -1) ) ; @#
            fclose(f); } @#
            
          @q  printf("%d  %d\n", rep, (X0 > 0 ? ngens : -1) ); @>

          rep = rep +  (X0 > 0.0 ? 1 : 0) ; }  @#
  
    

  @t \newline @>
  /* free memory */ 
    @t \newline @>

    free(Z); @#
      free(tXi); @#
      free(PXi); @#
        gsl_ran_discrete_free( Pmass); @#
   free( tempJuve); @#
   free( etimes); @#
   free(aindex); @#
  free( Pop); @#
  free( leffects); @#
 
}



@*1 {\bf run over parameters }. 


Run over some parameters. 

@<parameters@>=@#
void run_parameters(  int N, int nalleles,  double variance,     double a, double b,  int Psi, double
s, double znull, double epsilon, int nruns,   char skra[200],   gsl_rng * r  )
{
   double ai = 1. ; @#
   double out ; 
   int psii ; @#
   while( ai < 2. ){ @#
    for( psii = 100000; psii < 1000001 ; psii = psii + 100000 ){ @#
      out  =  new_simulator( N, nalleles, variance, (ai > 0. ? ai : 0.),  b,  psii, s , znull ,
      epsilon,  nruns, r);  @#

      printf( "%g %d %g\n", ai, psii, out); @#   

         FILE * f = fopen(skra, "a") ; @#
         
           fprintf(f,  "%g ", ai) ;  @#
           
           fprintf(f,  "%d ", psii) ;  @#
           fprintf(f,  "%g\n", out) ;  @#
            
            fclose(f); 
            
            } @#
   
     ai = ai + 0.1 ; } @#
}





 
@*1 {\bf the $main$ function}. 
\label{sec:main}

@C

@<Includes@>@#
@<random number generator@>@#
@<object definitions@>@#
@<initialize distribution for $X_i$@>@#
@<population update@>@#
@<replicates@>@#



int main(int argc, char * argv[])@#
{@#

 @t initialise the random number generator @>@#
 setup_rng( (unsigned long int)atoi(argv[1]) ) ; @#

@t \newline @>
/* the expected value $\mathbb{E}[X] =  11.09016 $ when  $(\alpha, \gamma) = (1.0, 
10^5)$; $\mathbb{E}[X] = 2.606051$ when  $(\alpha, \gamma) = (1.5, 10^5)$; 
$\mathbb{E}[X] = 1.644924$ when $(\alpha, \gamma) = (2.0, 10^5)$;  $\mathbb{E}[X] = 
 6.48 $ when $(\alpha, \gamma) = (1.0, 10^3)$;   */
@t \newline @>
@t \newline @>
 /* {\tt POP\_SiZE} is $N$; {\tt N\_ALLELES} is number of alleles; {\tt
 PSI\_TRUNCATION} is  $\gamma$ in  \eqref{eq:skew};  {\tt TRAIT\_OPTIMUM } is
 $z_0$ \eqref{eq:alg}, \eqref{eq:exp}; {\tt EPSILON} is $y$ the threshold frequency */
@t \newline @>
#define POP_SIZE 1000
#define N_ALLELES 100
#define VARIANCE 0.0
@t \newline @>
/*If {\tt ALPHA} $(\alpha)$ is 0 then the Poisson distribution is assumed with
mean  {\tt BETA} */
@t \newline @>
#define ALPHA 1.0
#define BETA 7.485471
#define PSI_TRUNCATION 100000
#define S_SELECTION 1.0
#define TRAIT_OPTIMUM 0.0
#define EPSILON 0.95
#define RUNS 10

  @t \newline @>
  simulator( POP_SIZE, N_ALLELES, VARIANCE, ALPHA, BETA, PSI_TRUNCATION, S_SELECTION,  TRAIT_OPTIMUM, EPSILON,  RUNS, argv[2],  rngtype);  @#
  @t \newline @>

 @q  run_parameters(  atoi(argv[1]),  atoi(argv[2]),  atof(argv[3]),      atof(argv[4]), atof(argv[5]), atoi(argv[6]),  atof(argv[7]) ,  atof(argv[8]), atof(argv[9]), atoi(argv[10]), argv[12],  rngtype );  @>
 
 @t \newline @>
gsl_rng_free( rngtype ) ; @#
 return GSL_SUCCESS ; @#
}



@* {\bf Includes}. 

@<Includes@>=@#
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_sort.h>
#include <assert.h>







@* {\bf References}. 


\bibliographystyle{genetics}
\bibliography{refs}







@* {\bf Funding}. 

Funded by DFG grant  325/17-1 to Wolfgang Stephan through   DFG  SPP Priority
Programme  1819: Rapid Evolutionary Adaptation (\url{https://dfg-spp1819.uni-hohenheim.de/en/105254}).  










@
\end{document}
