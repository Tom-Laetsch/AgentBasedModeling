#################################### GENERIC HELPER FUNCTIONS
safe_seq_len <- function( n ){
     tryCatch(
          { return( seq_len( n ) ) },
          error = function(e) return( integer(length = 0L) )
     )
}


dist1 <- function( dist, coef, min_val = 0 ){
     return( coef * max( dist, min_val ) )
}


dist2 <- function( dist, coef, min_val = 0 ){
     return( coef * max( dist, min_val )^2 )
}


dist3 <- function( dist, coef, min_val = 0 ){
     return( coef * max( dist, min_val )^3 )
}

distGraded <- function( dist, 
                        close_coef,
                        close_cutoff, 
                        mid_cutoff, 
                        min_val = 1 ){
     # coef chosen to create smooth cutoff boundaries
     mid_coef <- close_coef / close_cutoff
     far_coef <- mid_coef / mid_cutoff
     if( dist <= close_cutoff ){
          return( close_coef * dist )
     } else if( close_cutoff< dist & dist <= mid_cutoff ) {
          return( mid_coef * dist^2 )
     }
     return( far_coef * dist^3 )
}

SquareLabelledMatrixMaker <- function( labs ){
     n_dim <- length(labs)
     if( n_dim == 0 ) return( NULL )
     m <- matrix( rep(0, times = n_dim*n_dim), nrow = n_dim, ncol = n_dim )
     colnames( m ) <- labs
     rownames( m ) <- labs
     return( m )
}

PMFer_UniformByLabel <- function( labs ){
     n.labs <- length( labs )
     if( n.labs == 0 ) return( NULL )
     pmf <- rep(1, n.labs )
     names(pmf) <- labs 
     return( pmf/sum(pmf) )
}

PMFer_CombineWeightedPMFs <- function( wts, ... ){
     pmfs <- list( ... )
     n.pmfs <- length( pmfs )
     n.wts <- length( wts )
     if( n.pmfs == 0 ){
          warning("No PMFs found." )
          return( NULL )
     }
     if( n.wts == 0 ){
          warning("No weights found." )
          return( NULL )
     }
     if( n.wts != n.pmfs ){
          warning("Number of weights do not agree with with number of PMFs. Wrapping or truncation will occur.")
          while( n.pmfs > n.wts ){
               wts <- c( wts, wts )
               n.wts <- length( wts )
          }
          if( n.wts > n.pmfs ){
               wts <- wts[1:n.pmfs]
               n.wts <- n.pmfs
          }
     }
     if( any(wts < 0) ){
          warning("Some weights given are negative. Setting these to 0.")
          wts[wts < 0] <- 0
     }
     if( (den <- sum(wts)) == 0 ){
          warning( "Sum of the weights given equals 0. Setting all weights to 1." )
          wts[1:n.wts] <- 1
          den <- n.wts
     }
     wts_ps <- wts / den
     wt_pmf <- 0
     for( i in safe_seq_len(n.wts) ){
          if( i == 1 ){
               #used to make sure consistent naming
               names.pmf <- names(pmfs[[1]])
          }
          wt_pmf <- unname(wts_ps[i]) * pmfs[[i]][names.pmf] + wt_pmf
     }
     
     return( wt_pmf )
}

####################################### SECTION ###############################################
###############################################################################################
################################### HELPER FUNCTIONS ##########################################
###############################################################################################
###############################################################################################

RectNbhdMaker <- function( topl_r, topl_c, botr_r, botr_c ){
     "
     Given points topl_r, topl_c, botr_r, botr_c, creates
     a rectangular nbhd_layout object with top-left (topl_r, topl_c) 
     and bottom-right (botr_r, botr_c).
     "
     rlen <- botr_r - topl_r + 1
     clen <- botr_c - topl_c + 1
     bdr_n <- 2*(rlen+clen-2)
     intr_n <- rlen*clen - bdr_n
     
     bdr <- data.frame(node_type = rep('bdr', times = bdr_n))
     rows <- c( rep(topl_r, times = clen), rep(botr_r, times = clen) )
     if(rlen > 2) rows <- c(rows, rep( (topl_r+1):(botr_r-1), times = 2 ) )
     bdr$row <- as.integer( rows )
     cols <- rep(topl_c:botr_c, times = 2)
     cols <- c( cols, rep(topl_c, times = rlen-2), rep(botr_c, times = rlen-2) )
     bdr$col <- as.integer( cols )
     
     intr <- NULL
     if( intr_n > 0 ){
          intr <- data.frame(node_type = rep('intr', times = intr_n))
          intr$row <- as.integer( rep( (topl_r+1):(botr_r-1), each = clen-2 ) )
          intr$col <- as.integer( rep( (topl_c+1):(botr_c-1), times = rlen-2) )
     }
     
     return( rbind(bdr,intr) )
}

RectLatticeMaker <- function( nrows, ncols ){
     "
     Creates a rectangular lattice with specified num rows and cols
     "
     row <- rep( 1:nrows, each = ncols )
     col <- rep( 1:ncols, times = nrows )
     
     return( as.data.frame(cbind(row,col)) )
     
}


InBounds <- function( node, grid ){
     return( any(grid[,'row'] == node['row'] & grid[,'col'] == node['col']) )
}

rviolence <- function( n = 1, 
                       max_violence = 9, 
                       min_violence = 1, 
                       mean_violence = 3,
                       scale = 0.1 ){
     # Generate random violence factors for the agents
     # currently we generate from a binomial-like distribution
     N <- max_violence - min_violence
     binom_mean <- mean_violence - min_violence
     p <- binom_mean / N
     x <- rbinom( n = n, size = N, prob = p)
     
     return( scale * (x + min_violence ) )
}


HostAgentCntsHelper <- function( hostiles, qts = c(1.5, 1.25, 1, .75, .6) ){
     if( !exists('PLAYGROUND') ){
          stop( "PLAYGROUND must exist to use HostAgentHelper" )
     }
     if( !exists('.NBHDSUM') ){
          warning('.NBHDSUM not found. Creating .NBHDSUM..')
          if( !CreateNbhdSummary() ){
               stop( 'ERROR: Unable to creat NBHDSUM.')
          }
     }
     
     
     host_nbhds <- .NBHDSUM[.NBHDSUM$hostile, ]
     n.qts <- length(qts)
     if( n.qts < 1){
          warning( 'No quantile information passed. No changes made.' )
          return(hostiles)
     }
     if( n.qts >= 1 ){
          qs <- quantile( host_nbhds$n_nodes, probs = seq(0,1,1/n.qts) )
          n.qs <- length( qs )
     }
     
     get_fac <- function( n_nodes ){
          if( n_nodes <= qs[1] ) return( qts[1] )
          if( n_nodes >= qs[n.qs]) return( qts[n.qts] )
          for( i in 1:(n.qs-1) ){
               if( n_nodes > qs[i] & n_nodes <= qs[i+1] ){
                    return( qts[i] )
               }
          }
     }
     
     for( idx in safe_seq_len(nrow(host_nbhds)) ){
          info <- host_nbhds[idx, ]
          fac <- get_fac( info$n_nodes )
          hostiles[[info$nbhd_name]]$agents <- c(standard = floor(fac * info$n_nodes))
     }
     return( hostiles )
}


############################## EXISTENCE CHECKERS #############################################

PgExistChecker <- function( create_if_missing = T,
                            check_names = c('row',
                                            'col',
                                            'node_type',
                                            'designation',
                                            'nbhd_name',
                                            'family',
                                            'hostile')
                            ){
     pg_exist <- F
     if( exists('PLAYGROUND') ){
          if( !is.null(PLAYGROUND) ){
               if( 'data.frame' %in% class(PLAYGROUND) ){
                    if( all(check_names %in% names(PLAYGROUND)) ){
                         if( nrow(PLAYGROUND) > 0 ){
                              pg_exist <- T
                         }
                    }
               }
          }
     } else if( create_if_missing ){
          # try to create, but don't allow an infinite loop
          if( CreatePlayground() ) pg_exist <- PgExistChecker( F )
     }
     return( pg_exist )
}

AuthDFExistChecker <- function( create_if_missing = T, 
                                check_names = c('assigned_nbhd',
                                                'current_node',
                                                'uid')
                                ){
     auth_exist <- F
     if( exists('AUTHORITIES_DF') ){
          if( !is.null(AUTHORITIES_DF) ){
               if( 'data.frame' %in% class(AUTHORITIES_DF) ){
                    if( all(check_names %in% names(AUTHORITIES_DF)) ){
                         if( nrow(AUTHORITIES_DF) > 0 ){
                              auth_exist <- T
                         }
                    }
               }
          }
     } else if( create_if_missing ){
          # try to create, but don't allow an infinite loop
          if( CreateAuthoritiesDF() ) auth_exist <- AuthDFExistChecker( F )
     }
     return( auth_exist )
}

HostDFExistChecker <- function( create_if_missing = T,
                                check_names = c('designation',
                                                'family',
                                                'base_nbhd',
                                                'current_nbhd',
                                                'current_nbhd_desig',
                                                'base_node',
                                                'current_node',
                                                'agent_class',
                                                'violence_fac',
                                                'uid') 
                                ){
     host_exist <- F
     if( exists('HOSTILES_DF') ){
          if( !is.null(HOSTILES_DF) ){
               if( 'data.frame' %in% class(HOSTILES_DF) ){
                    if( all(check_names %in% names(HOSTILES_DF)) ){
                         if( nrow(HOSTILES_DF) > 0 ){
                              host_exist <- T
                         }
                    }
               }
          }
     } else if( create_if_missing ){
          # try to create, but don't allow an infinite loop
          if( CreateAuthoritiesDF() ) host_exist <- HostDFExistChecker( F )
     }
     return( host_exist )
}

############################### METRIC STUFF ##################################################

# distance between a node and a nbhd
NodeNbhdDist <- function( node, nbhd_coords ){
     return( min( apply( X=nbhd_coords, FUN=function(x){ METRIC( x, node ) }, MARGIN=1 ) ) )
}

# distance between two nbhds
NbhdNbhdDist <- function( nbhd_coords1, nbhd_coords2 ){
     return( min( apply( nbhd_coords1, MARGIN=1, FUN=function(x){ NodeNbhdDist( x, nbhd_coords2 ) } ) ) )
}

########################## NEIGHBOR CALCULATIONS #######################################


NodeToIndex <- function( node ){
     # requires PLAYGROUND
     w <- (PLAYGROUND$row == node['row']) & (PLAYGROUND$col == node['col'])
     return( rownames(PLAYGROUND[w, ]) )
}

IndexToNode <- function( idx ){
     # requires PLAYGROUND
     return( rownames(PLAYGROUND[idx, c('row','col')]) )
}

#tests if two nodes are neighbors 
AreNbrs <- function( node1, node2 ){
     # requires NBR_RADIUS
     return( METRIC(node1,node2) <= NBR_RADIUS )
}

NbrIndices <- function( node_index ){
     #requires PLAYGROUND
     
     if( exists('.NBR_INDICES_CACHE') ){
          if( !is.null(.NBR_INDICES_CACHE[[node_index]]) ){
               if( WALLED_NBHDS ){
                    return( .NBR_INDICES_CACHE[[node_index]]$walled )
               } else {
                    return( .NBR_INDICES_CACHE[[node_index]]$unwalled )
               }
          }
          nbr_indices_cache <- .NBR_INDICES_CACHE
          nbr_indices_cache[[node_index]] <- list()
     } else {
          nbr_indices_cache <- list()
          nbr_indices_cache[[node_index]] <- list()
     }
     
     node <- PLAYGROUND[node_index, ]
     row <- as.numeric(node['row'])
     col <- as.numeric(node['col'])
     diam <- 2*NBR_RADIUS + 1
     rows <- rep( (row-NBR_RADIUS):(row+NBR_RADIUS), each = diam )
     cols <- rep( (col-NBR_RADIUS):(col+NBR_RADIUS), times = diam )
     
     candidates <- data.frame( row = rows, col = cols )
     unwalled_nbr_idx <- unlist( apply( candidates, MARGIN = 1, 
                                        FUN = function( x ){
                                             ifelse( AreNbrs(x, node),
                                                     return( NodeToIndex( x ) ),
                                                     return( NULL ) )
                                        }) )
     nbr_indices_cache[[node_index]]$unwalled <- unwalled_nbr_idx
     pg_red <- PLAYGROUND[unwalled_nbr_idx, ]
     req_nbhd <- PLAYGROUND[node_index, ]$nbhd_name
     
     walled_nbr_idx <- rownames( pg_red[pg_red$nbhd_name == req_nbhd, ] )
     nbr_indices_cache[[node_index]]$walled <- walled_nbr_idx
     
     assign( '.NBR_INDICES_CACHE', nbr_indices_cache, .GlobalEnv )
     
     if( WALLED_NBHDS ){
          return( walled_nbr_idx )
     } else {
          return( unwalled_nbr_idx )
     }
     
}



################### Create* Family of Functions: Initialization ###############################

CreateNbhds <- function( ){
     # check that HOSTILES and OTHER_NBHDS are defined
     if( !( exists("HOSTILES") | exists("OTHER_NBHDS") ) ) return(F)
     
     # check that these are a list of lists
     combined_nbhds <- c(HOSTILES, OTHER_NBHDS)
     if( !all(unlist(lapply(combined_nbhds, FUN=function(x){class(x) == 'list'}))) ) return(F)
     
     nbhds <- list()
     # keep only wanted info and add the field "hostile" to indicate whehter a given nbhd is hostile
     host_nbhds <- names(HOSTILES)
     for( nbhd in names(combined_nbhds) ){
          cnbhd <- combined_nbhds[[nbhd]]
          nbhds[[nbhd]] <- list()
          nbhds[[nbhd]]$nbhd_lattice <- cnbhd$nbhd_lattice
          ifelse( is.null(cnbhd$designation),
                  nbhds[[nbhd]]$designation <- nbhd,
                  nbhds[[nbhd]]$designation <- cnbhd$designation )
          ifelse( is.null(cnbhd$family),
                  nbhds[[nbhd]]$family <- nbhd,
                  nbhds[[nbhd]]$family <- cnbhd$family )
          if( nbhd %in% host_nbhds ){
               nbhds[[nbhd]]$hostile <- T
          } else {
               nbhds[[nbhd]]$hostile <- F
          }
     }
     
     #store in our workspace as .NBHDS
     assign('.NBHDS', nbhds, .GlobalEnv)
     return( T )
}



CreatePlayground <- function( combine_adj = F ){
     
     if( exists("PLAYGROUND") ){
          ret <- readline( prompt = "PLAYGROUND already exists. Recreate? ([yes] or no): ")
          if( tolower(ret) == 'n' | tolower(ret) == 'no' ) return( T )
     }
     
     if( VERBOSE ) {
          time.start <- Sys.time()
          cat("Importing playground design...")
     }
     playground <- NULL
     
     # STEP 1: join the named nbhds
     nbhd_names <- names( .NBHDS )
     for( nbhd_name in nbhd_names ){
          nbhd <- .NBHDS[[nbhd_name]]
          nbhd_lattice <- nbhd$nbhd_lattice
          n.row <- nrow(nbhd_lattice)
          nbhd_lattice$designation <- nbhd$designation #rep( nbhd$designation, n.row )
          nbhd_lattice$nbhd_name <- nbhd_name #rep( nbhd_name, n.row )
          ifelse( is.null(nbhd$family),
                  nbhd_lattice$family <- nbhd_name,
                  nbhd_lattice$family <- nbhd$family
          )
          nbhd_lattice$hostile <- nbhd$hostile #rep( nbhd$hostile, n.row )
          playground <- rbind( playground, nbhd_lattice )
     }
     # end STEP 1
     
     # STEP 2: join the ether
     playground <- merge( playground, .LATTICE[,c('row','col')], by = c('row','col'), all.x = T, all.y = T)
     playground$designation[ is.na(playground$designation) ] <- "ether"
     playground$nbhd_name[ is.na(playground$nbhd_name) ] <- "ether"
     playground$family[ is.na(playground$family) ] <- "ether"
     playground$hostile[ is.na(playground$hostile) ] <- FALSE
     playground$node_type[ is.na(playground$node_type) ] <- "intr"
     # end STEP 2
     
     # STEP 3: give meaningful rownames
     rnames <- rownames( playground )
     for( area in unique( playground$nbhd_name ) ){
          l <- playground$nbhd_name == area
          rnames[ l ] <- paste0( area, "_node", 1:sum( l ) )
     }
     rownames( playground ) <- rnames
     # end STEP 3
     if( VERBOSE ) cat(" Done \n \n")
     
     # STEP 4: calculate inter-area distances 
     if( VERBOSE ) cat('Calculating relative distances. This can take time...\n')
     
     nbhds <- unique( playground$nbhd_name )
     non_ether_nbhds <- nbhds[ !('ether' == nbhds) ]
     n.nbhds <- length( nbhds )
     n.non_ether_nbhds <- length( non_ether_nbhds )
     n.row <- nrow( playground )
     
     ## initialize the "distance to" columns
     for( a in nbhds ){
          dist_from_a = paste0('dist_from_',a)
          playground[[dist_from_a]] <- numeric( length = n.row )
     }
     
     ## the following for loop is where the dist calc magic happens
     ether_nodes <- playground[playground$nbhd_name == 'ether', c('row','col')]
     for( an in safe_seq_len(n.non_ether_nbhds) ){
          if( VERBOSE ) inter.time.start <- Sys.time()
          a <- non_ether_nbhds[an]
          dist_from_a <- paste0('dist_from_',a)
          nbhd <- playground[playground$nbhd_name == a, c('row','col')]
          nbhd_bdr <- playground[playground$nbhd_name == a & playground$node_type == 'bdr', c('row','col')]
          #fill in distance to a from ether
          playground[[dist_from_a]][playground$nbhd_name == 'ether'] <- apply( ether_nodes, 
                                                                               MARGIN = 1, 
                                                                               FUN = function( b ){ 
                                                                                    NodeNbhdDist( b, nbhd_bdr ) 
                                                                               } 
          )
          ## ether is a special case since it surrounds each nbhd, so we can ...
          ### speed up calculations for distance to ether by only checking distance ...
          ### to bounding ether nodes (nearby_ether)
          ether <- playground[ playground$nbhd_name == 'ether', ]
          nearby_ether <- ether_nodes[ ether[dist_from_a] == 1, ]
          playground$dist_from_ether[playground$nbhd_name == a] <- apply( nbhd, 
                                                                          MARGIN = 1,
                                                                          FUN = function( b ){
                                                                               NodeNbhdDist( b, nearby_ether )
                                                                          }
          )
          ## if here an == n.non_ether_nbhds, we've iterated through everything and are done
          if( an == n.non_ether_nbhds ) break
          
          ## get relative distances between area a and other areas a2
          for( an2 in (an+1):n.non_ether_nbhds ){
               a2 <- non_ether_nbhds[an2]
               dist_from_a2 <- paste0('dist_from_',a2)
               nbhd2 <- playground[playground$nbhd_name == a2, c('row','col')]
               nbhd2_bdr <- playground[playground$nbhd_name == a2 & playground$node_type == 'bdr', c('row','col')]
               #fill in distance to a2 from a
               playground[[dist_from_a2]][playground$nbhd_name == a] <- apply( nbhd, 
                                                                               MARGIN = 1, 
                                                                               FUN = function( b ){ 
                                                                                    NodeNbhdDist( b, nbhd2_bdr ) 
                                                                               } 
               )
               ## fill in distance to a from a2
               playground[[dist_from_a]][playground$nbhd_name == a2] <- apply( nbhd2, 
                                                                               MARGIN = 1, 
                                                                               FUN = function( b ){ 
                                                                                    NodeNbhdDist( b, nbhd_bdr ) 
                                                                               } 
               )
          }
          if( VERBOSE ) {
               inter.time.end <- Sys.time()
               cat( sprintf("\t (%d/%d) finished with %s in %s\n",  an, n.non_ether_nbhds, a,
                            format(inter.time.end-inter.time.start)) )
          }
     }
     if( VERBOSE ) {
          inter.time.end <- Sys.time()
          cat( sprintf("\t (%d/%d) finished with %s in %s\n", an, n.non_ether_nbhds, a, 
                       format(inter.time.end-inter.time.start)) )
     }
     if( VERBOSE ) cat(" Done \n \n")
     # end STEP 4
     
     
     # STEP 5: Saving the fruit of our labours 
     if( VERBOSE ) cat("Saving the playground as PLAYGROUND in workspace...")
     assign('PLAYGROUND',playground,.GlobalEnv)
     if( VERBOSE ) cat(" Done \n \n")
     # end STEP 5
     
     # completion message
     if( VERBOSE ) {
          time.end <- Sys.time()
          cat( sprintf("Playground Finished. Total time: %s\n", format(time.end-time.start)) )
     }
     
     return( T )
}


CreatePlaygroundSummaries <- function( ){
     if( !PgExistChecker() ){
          warning('PLAYGROUND not found or empty. Cannot create .PGSUM.')
          return( F )
     }
     
     ## .PGSUM
     nb_grp <- dplyr::group_by(PLAYGROUND, designation, family, nbhd_name )
     pgsum <-   dplyr::summarise_at(nb_grp, .cols = dplyr::vars(dplyr::contains("dist_from_")), .funs = min )
     pgsum$hostile <- as.logical(
          ( dplyr::summarise(nb_grp, h = max(hostile)) )$h
     )
     assign('.PGSUM', as.data.frame( pgsum ), .GlobalEnv )
     
     ## .ETHSUM
     ether <- PLAYGROUND[PLAYGROUND$designation == "ether", ]
     host_desigs <- unique( PLAYGROUND$designation[PLAYGROUND$hostile] )
     ethsum <- NULL
     for( des in host_desigs ){
          dist_from_des <- paste0( "dist_from_", des )
          dist_from_cols <- paste0( "dist_from_", unique(PLAYGROUND$nbhd_name[PLAYGROUND$designation == des]) )
          temp_df <- as.data.frame( apply( ether[,dist_from_cols], MARGIN = 1, FUN = min ) )
          colnames( temp_df ) <- dist_from_des
          ifelse( is.null(ethsum), ethsum <- temp_df, ethsum <- cbind( ethsum, temp_df ) )
     }
     assign('.ETHSUM', ethsum, .GlobalEnv)
     
     return( T )
}

CreateHostilesDF <- function( ){
     
     #make sure everything is well defined
     if( !exists('HOSTILES') ){
          cat(paste("ERROR: HOSTILES list not found.","There is no default list for HOSTILES.\n\n", sep = " "))
          return( F )
     }
     if( !PgExistChecker() ){
          warning('PLAYGROUND not found or empty. Cannot create HOSTILES_DF')
          return( F )
     }
     
     host_nbhd_names <- names( HOSTILES )
     if( !all( host_nbhd_names %in% unique(PLAYGROUND$nbhd_name)) ){
          cat(paste("Warning: Some hostile nbhd names not represented in the playground.",
                    "Those not represented will not be assigned agents.\n\n", sep = " "))
     }
     if(VERBOSE) cat('Importing HOSTILES information... ')
     designation <- character()
     family <- character()
     agent_class <- character()
     current_node <- character()
     current_nbhd <- character()
     for( nbhd_name in host_nbhd_names ){
          nbhd <- HOSTILES[[nbhd_name]]
          rowns <- rownames( PLAYGROUND[PLAYGROUND$nbhd_name == nbhd_name, ] )
          agents <- nbhd$agents
          aclasses <- names(agents)
          for( aclass in aclasses ){
               n <- agents[[aclass]]
               ifelse( is.null( nbhd$designation ),
                       designation <- c( designation, rep(nbhd_name, times = n) ),
                       designation <- c( designation, rep(nbhd$designation, times = n) ) )
               ifelse( is.null( nbhd$family ),
                       family <- c( family, rep(nbhd_name, times = n) ),
                       family <- c( family, rep(nbhd$family, times = n) ) )
               agent_class <- c( agent_class, rep(aclass, times = n) )
               current_node <- c( current_node, sample( rowns, size = n, replace = TRUE ) )
               current_nbhd <- c( current_nbhd, rep(nbhd_name, times = n) )
          } 
     }
     hdf <- data.frame(   designation = designation,
                          family = family,
                          base_nbhd = current_nbhd,
                          current_nbhd = current_nbhd,
                          current_nbhd_desig = designation,
                          base_node = current_node,
                          current_node = current_node,
                          agent_class = agent_class,  
                          stringsAsFactors = FALSE )
     n.row <- nrow(hdf)
     hdf$violence_fac <- rviolence( n = n.row )
     # give a unique id to each hostile agent
     hdf$uid <- paste0( 'hostile_', 1:n.row )
     if(VERBOSE) cat('Done\n\n')
     
     if(VERBOSE) cat('Saving to HOSTILES_DF... ')
     # save HOSTILES_DF in our workspace
     assign('HOSTILES_DF', hdf, .GlobalEnv)
     if(VERBOSE) cat('Done\n\n')
     
     return( T )
}

CreateHostSummary <- function( ){
     if( !HostDFExistChecker() ){
          warning(paste("HOSTILES_DF not found. Cannot create .HOSTSUM."))
          return( F )
     }
     grp <- dplyr::group_by( HOSTILES_DF, designation, family, base_nbhd )
     df <- dplyr::summarise( grp, n_agents = n() )
     df <- as.data.frame( lapply( df, FUN = as.character ), stringsAsFactors = F )
     df$n_agents <- as.integer( df$n_agents )
     assign('.HOSTSUM', df, .GlobalEnv)
     return( T )
}

CreateNbhdSummary <- function( ){
     if( !PgExistChecker() ){
          warning('PLAYGROUND not found or empty. Cannot create .NBHDSUM.')
          return( F )
     }
     grp <- dplyr::group_by( PLAYGROUND, hostile, designation, family, nbhd_name )
     df <- dplyr::summarise( grp, n_nodes = n() ) 
     df <- as.data.frame( lapply( df, FUN = as.character ), stringsAsFactors = F )
     df$n_nodes <- as.integer( df$n_nodes )
     df$hostile <- as.logical( df$hostile )
     assign('.NBHDSUM', df, .GlobalEnv)
     return( T )
     
}

CreateAuthNbhdNodes <- function( ){
     "
     Creates a list, named by neighborhoods, indicating which nodes are within
     an authority agent's range of motion. This will be all nodes within the nbhd
     itself as well as all those within the ether of a certian radius specified
     within the AUTHORITIES list.
     "
     
     #make sure everything is well defined
     if( !PgExistChecker() ){
          warning('PLAYGROUND not found or empty. Cannot create AUTH_NBHD_NODES.')
          return( F )
     }
     if( !exists('AUTHORITIES') ){ 
          assign('AUTH_NBHD_NODES', list(), .GlobalEnv)
          return( T ) 
     }
     auth <- AUTHORITIES
     if( is.null(auth$assignment_radius) ) auth$assignment_radius <- 5
     nbhds <- PLAYGROUND[PLAYGROUND$nbhd_name != "ether", ]
     ether <- PLAYGROUND[PLAYGROUND$nbhd_name == "ether", ]
     auth_nbhd_nodes <- list()
     for( nbhd_name in unique(nbhds$nbhd_name) ){
          auth_nbhd_nodes[[nbhd_name]] <- rownames( nbhds[nbhds$nbhd_name == nbhd_name, ] )
          l <- ether[[paste0('dist_from_',nbhd_name)]] <= auth$assignment_radius
          auth_nbhd_nodes[[nbhd_name]] <- c( auth_nbhd_nodes[[nbhd_name]], rownames(ether[l, ]) )
     }
     
     assign('AUTH_NBHD_NODES', auth_nbhd_nodes, .GlobalEnv)
     
     return( T )
     
}


CreateEmptyAuthoritiesDF <- function(){
     #df <- data.frame( assigned_nbhd = character(), current_node = character() )
     assign('AUTHORITIES_DF', NULL, .GlobalEnv)
     return( T )
}

CreateAuthoritiesDF <- function( ){
     # Requires AUTHORITIES, AUTH_NBDH_NODES, PLAYGROUND, .PGSUM, .HOSTSUM
     
     if( !exists('AUTHORITIES') ){
          if( VERBOSE ){
               cat(paste("AUTHORITIES list not found.","Defaulting to no authority agents.\n\n", sep = " "))
          }
          CreateEmptyAuthoritiesDF()
          return( T )
     }
     
     auth <- AUTHORITIES
     n_free <- auth$n_free
     n_assigned <- auth$n_assigned
     
     if( is.null(n_free) | (n_free < 1) ) n_free <- 0
     if( is.null(n_assigned) | (n_assigned < 1) ) n_assigned <- 0
     
     if( n_free == 0 & n_assigned == 0 ){
          CreateEmptyAuthoritiesDF()
          return( T )
     }
     
     if( is.null(auth$min_per_nbhd) ) auth$min_per_nbhd <- 0
     if( is.null(auth$assignement_radius) ) auth$assignment_radius <- 0
     
     given_nbhds <- intersect( names(AUTHORITIES), .NBHDSUM$nbhd_name )
     remain_nbhds <- .NBHDSUM$nbhd_name[ !(.NBHDSUM$nbhd_name %in% given_nbhds) ]
     
     assigned_nbhds <- character()
     current_node <- character()
     if( length(given_nbhds) > 0 ){
          for( nbhd_name in given_nbhds ){
               if( !is.null( n <- auth[[nbhd_name]]$n ) ){
                    assigned_nbhds <- c( assigned_nbhds, rep(nbhd_name, times = n) )
                    current_node <- sample( AUTH_NBHD_NODES[[nbhd_name]], size = n, replace = T )
               }
          }
     }
     length.assigned <- length(assigned_nbhds)
     auth$n_assigned <- max( n_assigned, length.assigned )
     
     n.auth.agents.remain <- auth$n_assigned - length.assigned
     # remaining hostile nbhds 
     remain_nbhds <- .HOSTSUM$base_nbhd[ !(.HOSTSUM$base_nbhd %in% assigned_nbhds) ]
     n.nbhds.remain <- length( remain_nbhds )
     if( n.nbhds.remain == 0 ){
          auth$n_assigned <- length.assigned
          n_free <- auth$n_free + n.auth.agents.remain
     } else {
          n.auth.agents.reserved <- auth$min_per_nbhd*n.nbhds.remain
          n.auth.agents.remain <- max( n.auth.agents.remain,  n.auth.agents.reserved)
          wts <- .HOSTSUM$n_agents[ .HOSTSUM$base_nbhd %in% remain_nbhds ]
          if( (den <- sum(wts)) == 0 ){ 
               n.wts <- length(wts)
               wts <- rep( 1, times = n.wts )
               wts <- wts/n.wts
          } else {
               wts <- wts / den 
          }
          agent_partition <- sapply( wts, 
                                     FUN = function( x ){
                                          floor( x * (n.auth.agents.remain - n.auth.agents.reserved ) ) + auth$min_per_nbhd } 
          )
          names( agent_partition ) <- remain_nbhds
          # add whichever remain as free agents
          n_free <- auth$n_free +  n.auth.agents.remain - sum( agent_partition )
          for( nbhd_name in remain_nbhds ){
               assigned_nbhds <- c( assigned_nbhds, rep(nbhd_name, times = agent_partition[nbhd_name]) )
               current_node <- c( current_node, 
                                  sample( AUTH_NBHD_NODES[[nbhd_name]], size = agent_partition[nbhd_name], replace = T) )
          }
     }
     
     assigned_nbhds <- c( assigned_nbhds, rep("free_", times = n_free) )
     rem_nodes <- rownames( PLAYGROUND[!(PLAYGROUND$nbhd_name %in% assigned_nbhds), ] )
     if( length(rem_nodes) < n_free ){
          rem_nodes <- rownames(PLAYGROUND)
     }
     current_node <- c( current_node, sample( rem_nodes, size = n_free, replace = T) )
     
     df <- data.frame( assigned_nbhd = assigned_nbhds, current_node = current_node, 
                       stringsAsFactors = F )
     
     # if there are any authorities, give them uids     
     if((n.df <- nrow(df)) > 0) df$uid <- paste0('auth_',safe_seq_len(n.df))
     
     assign('AUTHORITIES_DF', df, .GlobalEnv)
     assign('AUTHORITIES', auth, .GlobalEnv)
     
     return( T )
}

CreateAuthSummary <- function( ){
     if( !AuthDFExistChecker() ){
          warning('AUTHORITIES_DF not found; .AUTHSUM set to NULL.')
          assign('.AUTHSUM',NULL,.GlobalEnv)
          return(F)
     } 
     grp <- dplyr::group_by( AUTHORITIES_DF, assigned_nbhd )
     df <- as.data.frame( dplyr::summarise( grp, n_agents = n() ) )
     assign('.AUTHSUM', df, .GlobalEnv)
     return( T )
}

CreateDominanceMatrix <- function( ){
     # requires HOSTILES_DF
     dm <- SquareLabelledMatrixMaker( labs = unique( HOSTILES_DF$family ) )
     if( is.null(dm) ){
          cat('Failed to create Dominance Matrix.\n')
          return(F)
     }
     assign('DOMINANCE_MATRIX', dm, .GlobalEnv)
     return(T)
}


CreateTensionMatrix <- function( ){
     if( !HostDFExistChecker() ){
          warning('Must create HOSTILES_DF before TENSION_MATRIX')
          return(F)
     }
     if( !exists('MIN_TENSION') ){
          warning('MIN_TENSION value not found. Assigning default 0.2')
          assign('MIN_TENSION', 0.2, .GlobalEnv)
     } 
     labs <- unique( HOSTILES_DF$base_nbhd )
     tm <- SquareLabelledMatrixMaker( labs = labs )
     if( is.null(tm) ){
          cat('Failed to create Tension Matrix.\n')
          return(F)
     }
     
     for( nbhd1 in labs ){
          for( nbhd2 in labs ){
               if( NbhdDesignation(nbhd1) != NbhdDesignation(nbhd2) ){
                    tm[nbhd1,nbhd2] <- MIN_TENSION
               }
          }
     }
     
     assign('TENSION_MATRIX', tm, .GlobalEnv)
     return(T)
}


############################### SUB-SECTION: INTERACTION BEHAVIOR ###################
#####################################################################################

###################################### Update* Family of Functions

UpdateDominanceMatrix <- function( timestep ){
     if( !exists('DOMINANCE_MATRIX') ){
          if( !CreateDominanceMatrix() ){
               stop('Unable to create DOMINANCE_MATRIX. Cannot continue.')
          }
     }
     if(!exists('HOSTILITY_HIST')) return( F )
     dm <- DOMINANCE_MATRIX
     rnames <- rownames(dm)
     cnames <- colnames(dm)
     
     cur_hist <- HOSTILITY_HIST[HOSTILITY_HIST$timestep == timestep,]
     if( nrow(cur_hist) < 1 ) return(T)
     
     cnts <- as.data.frame( dplyr::summarise( dplyr::group_by(cur_hist, perp_fam, vict_fam ), n = n()) )
     cnts <- merge( x = cnts, 
                    y = data.frame(perp_fam = cnts$vict_fam, vict_fam = cnts$perp_fam ), 
                    by = c('perp_fam','vict_fam'),
                    all.x = T,
                    all.y = T )
     cnts$n[ is.na( cnts$n ) ] <- 0
     
     for( famn in safe_seq_len(nrow(cnts)) ){
          fi <- cnts[famn, ]$perp_fam
          fj <- cnts[famn, ]$vict_fam
          i_to_j <- cnts[ (cnts$perp_fam == fi) & (cnts$vict_fam == fj),  ]$n
          j_to_i <- cnts[ (cnts$perp_fam == fj) & (cnts$vict_fam == fi),  ]$n
          dm[fi == rnames, fj == cnames] <- DOM_MATRIX_ENTRY( outward = i_to_j, 
                                                              inward = j_to_i, 
                                                              old_dom = dm[fi == rnames, fj == cnames] )
     }
     assign('DOMINANCE_MATRIX', dm, .GlobalEnv)
     
     return(T)
}

UpdateTensionMatrix <- function( timestep ){
     if( !exists('TENSION_MATRIX') ){
          if( !CreateTensionMatrix() ){
               stop('Unable to create TENSION_MATRIX. Cannot continue.')
          }
     }
     tm <- TENSION_MATRIX
     rnames <- rownames(tm)
     cnames <- colnames(tm)
     
     if(!exists('HOSTILITY_HIST')) return( F )
     cur_hist <- HOSTILITY_HIST[HOSTILITY_HIST$timestep == timestep,]
     if( nrow(cur_hist) < 1 ) return(T)
     
     no_atk_cnts <- cur_hist[!(cur_hist$during_attack), ]
     if( nrow(no_atk_cnts) >  0 ){
          cnts <- as.data.frame( dplyr::summarise( dplyr::group_by(no_atk_cnts, 
                                                                   perp_base_nbhd, 
                                                                   vict_base_nbhd ), 
                                                   n = n()),
                                 stringsAsFactors = F)
     } else {
          cnts <- data.frame()
     }
     
     atk_cnts <- cur_hist[cur_hist$during_attack, ]
     redo_cnts_grp <- F
     if( nrow(atk_cnts) > 0 ){
          redo_cnts_grp <- T
          atk_cnts <- as.data.frame( dplyr::summarise( dplyr::group_by( atk_cnts, 
                                                                        perp_fam, 
                                                                        vict_fam ), 
                                                       n = n()),
                                     stringsAsFactors = F )
          
          
          #append to cnts so that effectively each nbhd involved in the attack feels the effect
          for( i in safe_seq_len(nrow(atk_cnts)) ){
               atk <- atk_cnts[i, ]
               perp_nbhds <- FamNbhd( atk$perp_fam )
               vict_nbhds <- FamNbhd( atk$vict_fam )
               for( pn in perp_nbhds ){
                    for( vn in vict_nbhds ){
                         cnts <- rbind( cnts, 
                                        data.frame( perp_base_nbhd = pn, vict_base_nbhd = vn, n = atk$n ),
                                        make.row.names = F,
                                        stringsAsFactors = F )
                    }
               }
          }
     }
     
     if( redo_cnts_grp ) {
          cnts <- as.data.frame( dplyr::summarise( dplyr::group_by(cnts, 
                                                                   perp_base_nbhd, 
                                                                   vict_base_nbhd ), 
                                                   n = sum(n)),
                                 stringsAsFactors = F)
     }
     
     # symmetrize the perp/vict nbhds for ease of calculation below
     cnts <- merge( x = cnts, 
                    y = data.frame(perp_base_nbhd = cnts$vict_base_nbhd, vict_base_nbhd = cnts$perp_base_nbhd ), 
                    by = c('perp_base_nbhd','vict_base_nbhd'),
                    all.x = T,
                    all.y = T )
     cnts$n[ is.na( cnts$n ) ] <- 0
     
     for( nbhdn in safe_seq_len(nrow(cnts)) ){
          ni <- cnts[nbhdn, ]$perp_base_nbhd
          nj <- cnts[nbhdn, ]$vict_base_nbhd
          i_to_j <- cnts[ (cnts$perp_base_nbhd == ni) & (cnts$vict_base_nbhd == nj),  ]$n
          j_to_i <- cnts[ (cnts$perp_base_nbhd == nj) & (cnts$vict_base_nbhd == ni),  ]$n
          
          tm[ni == rnames, nj == cnames] <- TENSION_MATRIX_ENTRY( outward = i_to_j, 
                                                                  inward = j_to_i, 
                                                                  old_tens = tm[ni == rnames, nj == cnames] )
          
     }
     assign('TENSION_MATRIX', tm, .GlobalEnv)
     
     return(T)
}

UpdateAttackHostilities <- function( family, 
                                     to_attack, 
                                     timestep ){
     
     if( !exists('DOMINANCE_MATRIX') ){
          if( !CreateDominanceMatrix() ){
               stop( 'ERROR: DOMINANCE_MATRIX cannot be found nor created to update during attack.' )
          }
     }
     
     host_hist <- data.frame()
     for( fam_to_atk in to_attack ){
          atk_dom <- ATTACK_DOM_COUNT( DOMINANCE_MATRIX[family, fam_to_atk] )
          perp_fam <- rep( family, times = atk_dom )
          perp_desig <- FamDesignation( perp_fam )
          vict_fam <- rep( fam_to_atk, times = atk_dom )
          vict_desig <- FamDesignation( vict_fam )
          
          host_hist <- rbind( host_hist,
                              data.frame( perp_id = NA, 
                                          perp_in_node = NA,
                                          perp_base_nbhd = NA,
                                          perp_fam = perp_fam,
                                          perp_desig = perp_desig, 
                                          vict_id = NA,
                                          vict_in_node = NA,
                                          vict_base_nbhd = NA,
                                          vict_fam = vict_fam,
                                          vict_desig = vict_desig,
                                          during_attack = T, 
                                          timestep = timestep),
                              make.row.names = F,
                              stringsAsFactors = F )
     }
     
     if( exists('HOSTILITY_HIST') ) host_hist <- rbind( HOSTILITY_HIST, host_hist, make.row.names = F, stringsAsFactors = F )
     assign( 'HOSTILITY_HIST', host_hist, .GlobalEnv )
     
     return(T)
}

UpdateInteractions <- function( timestep ){
     
     if( !HostDFExistChecker() ){
          warning( "HOSTILES_DF not found or empty. No interactions calcuated." )
          return( NULL )
     }
     
     ## initialize the vectors that we want to remember
     ### for ENEMY_SIGHTING_HIST
     witness_id <- character()
     witness_in_node <- character()
     witness_base_nbhd <- character()
     witness_in_nbhd <- character()
     witness_fam <- character()
     witness_desig <- character()
     observed_id <- character()
     observed_in_node <- character()
     observed_base_nbhd <- character()
     observed_fam <- character()
     observed_desig <- character()
     ### for HOSTILITY_HIST
     perp_id <- character()
     perp_in_node <- character()
     perp_base_nbhd <- character()
     perp_fam <- character()
     perp_desig <- character()
     vict_id <- character()
     vict_in_node <- character()
     vict_base_nbhd <- character()
     vict_fam <- character()
     vict_desig <- character()
     during_attack <- logical()
     
     who_moved <- HOSTILES_DF$current_nbhd_desig != HOSTILES_DF$designation
     moved_agents <- HOSTILES_DF[ who_moved, ]
     unmoved_agents <- HOSTILES_DF[ !who_moved, ]
     
     n.moved <- nrow( moved_agents )
     for( mv_n in safe_seq_len(n.moved) ){
          agent <- moved_agents[mv_n, ]
          agent_nbrs <- as.character( agent$current_node )
          who_enemies <- moved_agents$current_node %in% agent_nbrs
          who_enemies <- who_enemies & (moved_agents$designation != agent$designation)
          moved_enemies <- moved_agents[ who_enemies, ] 
          n.moved_enemies <- nrow( moved_enemies )
          for( mv_n2 in safe_seq_len(n.moved_enemies) ){
               agent2 <- moved_enemies[mv_n2, ]
               # update the hostility encounter log
               witness_id <- c( witness_id, agent$uid)
               witness_in_node <- c( witness_in_node, agent$current_node )
               witness_base_nbhd <- c( witness_base_nbhd, agent$base_nbhd )
               witness_in_nbhd <- c( witness_in_nbhd, agent$current_nbhd )
               witness_fam <- c( witness_fam, agent$family )
               witness_desig <- c( witness_desig, agent$designation )
               observed_id <- c( observed_id, agent2$uid )
               observed_in_node <- c( observed_in_node, agent2$current_node )
               observed_base_nbhd <- c( observed_base_nbhd, agent2$base_nbhd )
               observed_fam <- c( observed_fam, agent2$family )
               observed_desig <- c( observed_desig, agent2$designation )
               if( SHOOT(perp_agent = agent, vict_agent = agent2) ){
                    # update the hostility log
                    perp_id <- c( perp_id, agent$uid )
                    perp_in_node <- c( perp_in_node, agent$current_node )
                    perp_base_nbhd <- c( perp_base_nbhd, agent$base_nbhd )
                    perp_fam <- c( perp_fam, agent$family )
                    perp_desig <- c( perp_desig, agent$designation )
                    vict_id <- c( vict_id, agent2$uid )
                    vict_in_node <- c( vict_in_node, agent2$current_node )
                    vict_base_nbhd <- c( vict_base_nbhd, agent2$base_nbhd )
                    vict_fam <- c( vict_fam, agent2$family )
                    vict_desig <- c( vict_desig, agent2$designation )
               }
          }
          
          who_enemies <- unmoved_agents$current_node %in% agent_nbrs
          who_enemies <- who_enemies & (unmoved_agents$designation != agent$designation)
          unmoved_enemies <- unmoved_agents[ who_enemies, ] 
          n.unmoved_enemies <- nrow( unmoved_enemies )
          for( mv_n2 in safe_seq_len(n.unmoved_enemies) ){
               agent2 <- unmoved_enemies[mv_n2, ]
               # update the hostility encounter log agent -> agent2
               witness_id <- c( witness_id, agent$uid)
               witness_in_node <- c( witness_in_node, agent$current_node )
               witness_base_nbhd <- c( witness_base_nbhd, agent$base_nbhd )
               witness_in_nbhd <- c( witness_in_nbhd, agent$current_nbhd )
               witness_fam <- c( witness_fam, agent$family )
               witness_desig <- c( witness_desig, agent$designation )
               observed_id <- c( observed_id, agent2$uid )
               observed_in_node <- c( observed_in_node, agent2$current_node )
               observed_base_nbhd <- c( observed_base_nbhd, agent2$base_nbhd )
               observed_fam <- c( observed_fam, agent2$family )
               observed_desig <- c( observed_desig, agent2$designation )
               # update the hostility encounter log agent2 -> agent
               witness_id <- c( witness_id, agent2$uid)
               witness_in_node <- c( witness_in_node, agent2$current_node )
               witness_base_nbhd <- c( witness_base_nbhd, agent2$base_nbhd )
               witness_in_nbhd <- c( witness_in_nbhd, agent2$current_nbhd )
               witness_fam <- c( witness_fam, agent2$family )
               witness_desig <- c( witness_desig, agent2$designation )
               observed_id <- c( observed_id, agent$uid )
               observed_in_node <- c( observed_in_node, agent$current_node )
               observed_base_nbhd <- c( observed_base_nbhd, agent$base_nbhd )
               observed_fam <- c( observed_fam, agent$family )
               observed_desig <- c( observed_desig, agent$designation )
               if( SHOOT(perp_agent = agent, vict_agent = agent2) ){
                    # update the hostility log agent -> agent2
                    perp_id <- c( perp_id, agent$uid )
                    perp_in_node <- c( perp_in_node, agent$current_node )
                    perp_base_nbhd <- c( perp_base_nbhd, agent$base_nbhd )
                    perp_fam <- c( perp_fam, agent$family )
                    perp_desig <- c( perp_desig, agent$designation )
                    vict_id <- c( vict_id, agent2$uid )
                    vict_in_node <- c( vict_in_node, agent2$current_node )
                    vict_base_nbhd <- c( vict_base_nbhd, agent2$base_nbhd )
                    vict_fam <- c( vict_fam, agent2$family )
                    vict_desig <- c( vict_desig, agent2$designation )
               }
               if( SHOOT(perp_agent = agent2, vict_agent = agent) ){
                    # update the hostility log agent2 -> agent
                    perp_id <- c( perp_id, agent2$uid )
                    perp_in_node <- c( perp_in_node, agent2$current_node )
                    perp_base_nbhd <- c( perp_base_nbhd, agent2$base_nbhd )
                    perp_fam <- c( perp_fam, agent2$family )
                    perp_desig <- c( perp_desig, agent2$designation )
                    vict_id <- c( vict_id, agent$uid )
                    vict_in_node <- c( vict_in_node, agent$current_node )
                    vict_base_nbhd <- c( vict_base_nbhd, agent$base_nbhd )
                    vict_fam <- c( vict_fam, agent$family )
                    vict_desig <- c( vict_desig, agent$designation )
               }
          }
     }
     
     host_sighting <- data.frame( witness_id = witness_id, 
                                  witness_in_node = witness_in_node,
                                  witness_base_nbhd = witness_base_nbhd,
                                  witness_in_nbhd = witness_in_nbhd,
                                  witness_fam = witness_fam,
                                  witness_desig = witness_desig,
                                  observed_id = observed_id,
                                  observed_in_node = observed_in_node,
                                  observed_base_nbhd = observed_base_nbhd,
                                  observed_fam = observed_fam,
                                  observed_desig = observed_desig,
                                  stringsAsFactors = F )
     
     
     host_hist <- data.frame( perp_id = perp_id,
                              perp_in_node = perp_in_node,
                              perp_base_nbhd = perp_base_nbhd,
                              perp_fam = perp_fam, 
                              perp_desig = perp_desig,
                              vict_id = vict_id,
                              vict_in_node = vict_in_node,
                              vict_base_nbhd = vict_base_nbhd,
                              vict_fam = vict_fam,
                              vict_desig = vict_desig,
                              stringsAsFactors = F  )
     
     if( nrow( host_hist ) > 0 ){
          # combine events where perps from the same nbhd, standing at the same node shoot at victs from the same nbhd
          host_hist <- dplyr::filter( .data = dplyr::group_by( .data = host_hist, perp_in_node, perp_base_nbhd, vict_base_nbhd ),
                                      row_number() == sample(safe_seq_len(n()), size = 1) )
          # combind events where multiple perps from the same nbhd shoot at the victims in the same node from the same nbhd
          host_hist <- dplyr::filter( .data = dplyr::group_by( .data = host_hist, vict_in_node, vict_base_nbhd, perp_base_nbhd ),
                                      row_number() == sample(safe_seq_len(n()), size = 1) )
          host_hist <- as.data.frame( lapply( host_hist, FUN = as.character ), stringsAsFactors = F )
     }
     if( nrow(host_hist) > 0 ){
          host_hist$timestep <- as.integer(timestep)
          host_hist$during_attack <- F
     } else {
          host_hist$timestep <- integer()
          host_hist$during_attack <- logical()
     }
     if( nrow( host_sighting ) > 0) {
          # combine sightings of agents from the same nbhd at the same node witnessing enemies from the same family
          host_sighting <- dplyr::filter( .data = dplyr::group_by( .data = host_sighting, witness_in_node, witness_base_nbhd, observed_fam ),
                                          row_number() == sample(safe_seq_len(n()), size = 1) )
          host_sighting <- as.data.frame( lapply( host_sighting, FUN = as.character ), stringsAsFactors = F  )
     }
     ifelse( nrow(host_sighting) > 0,
             host_sighting$timestep <- as.integer(timestep),
             host_sighting$timestep <- integer() )
     
     auth_sighting <- data.frame()
     if( AuthDFExistChecker() ){
          if( ALL_AUTH_INTERACTIONS ){
               auth_nbrs <- lapply(AUTHORITIES_DF$current_node, FUN=NbrIndices)
               if( (n.auth <- nrow(AUTHORITIES_DF)) == length(auth_nbrs) ){
                    for( n in safe_seq_len(n.auth) ){
                         nbrs <- auth_nbrs[[n]]
                         hosts <- subset( HOSTILES_DF, current_node %in% nbrs )
                         n.hosts <- nrow(hosts)
                         auth_df <- AUTHORITIES_DF[rep(n, times = n.hosts), ]
                         witness_id <- hosts$uid
                         witness_in_node <- hosts$current_node
                         witness_base_nbhd <- hosts$base_nbhd
                         witness_fam <- hosts$family
                         witness_desig <- hosts$desig
                         auth_sighting <- rbind( auth_sighting,
                                                 dplyr::mutate(auth_df, 
                                                               witness_id = witness_id,
                                                               witness_in_node = witness_in_node, 
                                                               witness_base_nbhd = witness_base_nbhd,
                                                               witness_fam = witness_fam,
                                                               witness_desig = witness_desig,
                                                               timestep = timestep),
                                                 make.row.names = F )
                    }
                    if( nrow(auth_sighting) > 0 ){
                         # check if the auth_sighting prevented a hostility
                         auth_sighting <- as.data.frame( lapply( auth_sighting, FUN = as.character ), stringsAsFactors = F )
                         auth_sighting$timestep <- as.integer( auth_sighting$timestep )
                         
                         auth_to_host <- auth_sighting$witness_id %in% host_hist$perp_id
                         if( AGGL_AUTHORITY ){
                              auth_to_host <- auth_to_host | (auth_sighting$witness_id %in% host_hist$vict_id)
                         }
                         auth_sighting$prevented_host <- auth_to_host
                         
                         host_to_auth <- host_hist$perp_id %in% auth_sighting$witness_id
                         if( AGGL_AUTHORITY ){
                              host_to_auth <- host_to_auth | (host_hist$vict_id %in% auth_sighting$witness_id)
                         }
                         host_hist <- host_hist[ !host_to_auth, ]
                    }
               }
          } else {
               n.host_hist <- nrow(host_hist)
               if( n.host_hist > 0 ){
                    for( hostn in safe_seq_len(n.host_hist) ){
                         host <- host_hist[hostn, ]
                         nbrs <- NbrIndices(as.character(host$perp_in_node))
                         if( AGGL_AUTHORITY ){
                              nbrs <- c( nbrs, NbrIndices(as.character(host$vict_in_node)) )
                         }
                         nearby_auth <- ( AUTHORITIES_DF$current_node %in% nbrs )
                         auth_df <- AUTHORITIES_DF[nearby_auth, ]
                         if( nrow(auth_df) > 0 ){
                              auth_sighting <- rbind(auth_sighting, 
                                                     dplyr::mutate( auth_df, 
                                                                    witness_id = as.character(host$perp_id),
                                                                    witness_in_node = as.character(host$perp_in_node),
                                                                    witness_base_nbhd = as.character(host$perp_base_nbhd),
                                                                    witness_fam = as.character(host$perp_fam),
                                                                    witness_desig = as.character(host$perp_desig),
                                                                    timestep = timestep,
                                                                    prevented_host = T ),
                                                     make.row.names = F )
                         }
                    }
                    if( nrow(auth_sighting) > 0 ){
                         # remove those prevented hostilities from the record (they never happened)
                         auth_sighting <- as.data.frame( lapply( auth_sighting, FUN = as.character ), stringsAsFactors = F )
                         auth_sighting$timestep <- as.integer( auth_sighting$timestep )
                         
                         host_to_auth <- host_hist$perp_id %in% auth_sighting$witness_id
                         host_hist <- host_hist[ !host_to_auth, ]
                    }
               }
          }
     }
     
     # if( nrow( auth_sighting ) > 0 ){
     #      # combine sightings of agents from the same nbhd at the same node witnessing authorities
     #      auth_sighting <- dplyr::filter( .data = dplyr::group_by( .data = auth_sighting, witness_in_node, witness_base_nbhd ),
     #                                      row_number() == sample(safe_seq_len(n()), size = 1) )
     # }
     
     
     if( nrow(host_hist) > 0 ){
          if( exists('HOSTILITY_HIST') ){
               host_hist <- rbind( HOSTILITY_HIST, host_hist, make.row.names = F )
          }
          assign( 'HOSTILITY_HIST', host_hist, .GlobalEnv )
     }
     if( nrow(host_sighting) > 0 ){
          if( exists('ENEMY_SIGHTING_HIST') ){
               host_sighting <- rbind( ENEMY_SIGHTING_HIST, host_sighting, make.row.names = F )
          }
          assign( 'ENEMY_SIGHTING_HIST', host_sighting, .GlobalEnv )
     }
     if( nrow( auth_sighting ) > 0 ){
          if( exists('AUTH_SIGHTING_HIST') ){
               auth_sighting <- rbind( AUTH_SIGHTING_HIST, auth_sighting, make.row.names = F )
          }
          assign( 'AUTH_SIGHTING_HIST', auth_sighting, .GlobalEnv )
     }
     
     return( T )
}

#################################### DATA IO ##########################################

SavePlayground <- function( file_path = NULL ){
     tryCatch( {
          saveRDS( PLAYGROUND, 
                   file = file_path )
     },
     error = function(e){
          file_path <- file.choose()
          if( file.exists(file_path) ){
               ret <- readline(sprintf('Overwrite file: %s? (yes or [no]): ',
                                       basename(file_path)))
               if( tolower(ret) == 'yes' | tolower(ret) == 'y' ){
                    SavePlayground( file_path )
               } else {
                    SavePlayground( NULL )
               } 
          } else {
               SavePlayground(file_path)
          }
          return(T)
     } 
     )
}

OpenPlayground <- function( file_path = NULL ){
     tryCatch( {
          pg <- readRDS(file_path)
          assign('PLAYGROUND', pg, .GlobalEnv)
          return(T)
     },
     error = function(e){
          file_path <- file.choose()
          return( OpenPlayground( file_path ) )
     } 
     )
}


#################################### VISUALIZATION ####################################
.PalettePre <- function( agent_df ){
     colrs <- c('goldenrod','steelblue','springgreen','tan','mediumpurple','olivedrab')
     desigs <- unique(agent_df$designation)
     n.desigs <- length( desigs )
     while( n.desigs > length(colrs) ) colrs <- c(colrs, colrs )
     ret_pal <- NULL
     for( i in safe_seq_len(n.desigs) ){
          desig <- desigs[i]
          df <- subset( agent_df, designation == desig )
          fams <- unique( df$family )
          colr <- colrs[i]
          getPalette = grDevices::colorRampPalette( c(colr, paste0(colr,1:4) ) )
          pal <- getPalette( length( fams ) )
          names(pal) <- fams
          ret_pal <- c(ret_pal, pal)
     }
     
     return(ret_pal)
}

.ViewAgentsPre <- function( agent_df ){
     if( !PgExistChecker() ) stop('PLAYGROUND not found or empty. Cannot create plot.')
     max.row <- max(PLAYGROUND$row)
     agent_plot <- agent_df
     agent_plot$row <- (max.row + 1 - PLAYGROUND[agent_plot$current_node, ]$row)
     agent_plot$col <- PLAYGROUND[agent_plot$current_node, ]$col
     return( agent_plot )
}

ViewPlayground <- function() {
     cat('ViewPlayground has been taken down until reworked. However, ViewHostiles, ViewAuthorities, and ViewAgents are up and running.\n')
}

ViewHostiles <- function( alpha = .75, only_spillover = F ){
     if( !PgExistChecker() ){ 
          warning('PLAYGROUND not found. Cannot plot agents.')
          return( ggplot2::ggplot() )
     }
     if( !HostDFExistChecker() ){ 
          warning('HOSTILES_DF not found. Cannot plot agents.')
          return( ggplot2::ggplot() )
     }
     hdf <- HOSTILES_DF
     if( only_spillover ){
          hdf <- subset( hdf, current_nbhd != base_nbhd )
     }
     hplot <- .ViewAgentsPre( hdf )
     
     pal <- .PalettePre( hdf )
     if( is.null(pal) ) return(NULL)
     
     g <- ggplot2::ggplot( hplot, ggplot2::aes(x = col, y = row, color = family, shape = designation) )
     g <- g + ggplot2::scale_color_manual( values = pal )
     g <- g + ggplot2::geom_point( alpha = alpha )
     return( g )
}

ViewAuthorities <- function( alpha = 1 ){
     if( !PgExistChecker() ) stop( 'PLAYGROUND not found. Cannot create plot.' )
     if( !AuthDFExistChecker() ){ 
          warning('AUTHORITIES_DF not found. Cannot plot agents.')
          return( ggplot2::ggplot() )
     }
     adf <- AUTHORITIES_DF
     aplot <- .ViewAgentsPre( adf )
     
     g <- ggplot2::ggplot( aplot, ggplot2::aes(x = col, y = row) )
     g <- g + ggplot2::geom_point( color = 'firebrick3', alpha = alpha )
     return( g )
}

ViewAgents <- function( halpha = .5, 
                        aalpha = 1, 
                        hsize = 1.5,
                        asize = .6,
                        only_spillover = F ){
     if( !PgExistChecker() ) stop( 'PLAYGROUND not found or empty. Cannot create plot.' )
     if( !AuthDFExistChecker() ){ 
          # only plotting hostile agents since AUTHROTIES_DF not found
          return( ViewHostiles( alpha = halpha, only_spillover = only_spillover ) )
     }
     if( !HostDFExistChecker() ){
          # AUTHORITIES_DF not found. Only plotting hostile agents
          return( ViewAuthorities( alpha = aalpha ) )
     }
     hdf <- HOSTILES_DF
     if( only_spillover ){
          hdf <- subset( hdf, current_nbhd != base_nbhd )
     }
     hplot <- .ViewAgentsPre( hdf )
     hplot$alpha <- 'hostile'
     hplot$size <- 'hostile'
     hplot <- dplyr::select(hplot, col, row, family, designation, alpha, size )
     
     adf <- AUTHORITIES_DF
     aplot <- .ViewAgentsPre( adf )
     aplot$designation <- 'authority'
     aplot$family <- 'authority'
     aplot$alpha <- 'authority'
     aplot$size <- 'authority'
     aplot <- dplyr::select(aplot, col, row, family, designation, alpha, size)
     
     p <- rbind( hplot, aplot )
     
     
     pal <- .PalettePre( hdf )
     if( is.null(pal) ) return(NULL)
     pal <- c(pal, authority = 'firebrick3')
     
     g <- ggplot2::ggplot( p, ggplot2::aes(x = col, y = row, color = family, shape = designation, alpha = alpha, size = size ) )
     g <- g + ggplot2::scale_color_manual( values = pal )
     g <- g + ggplot2::scale_alpha_manual( values = c(hostile = halpha, authority = aalpha), guide = F )
     g <- g + ggplot2::scale_size_manual( values = c(hostile = hsize, authority = asize), guide = F )
     g <- g + ggplot2::geom_point(  )
     return( g )
}

.VisPlayground <- function( alph = 1 ){
     if( !exists('PLAYGROUND') ){
          cat('Must create PLAYGROUND using CreatePlayground before visualization.\n')
          return( F )
     }
     if( !all( c('ggplot2','dplyr','RColorBrewer') %in% installed.packages() ) ){
          cat('Necessary packages for ViewPlayground are ggplot2, RColorBrewer, and dplyr.\n')
          return(F)
     }
     
     max_row <- max(PLAYGROUND$row)
     pgplot <- PLAYGROUND[,c('row','col','designation' )]
     pgplot <- dplyr::mutate(pgplot, plot_row = max_row + 1 - row )
     
     colourCount = length(unique(pgplot$designation))
     getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))
     
     g <- ggplot2::ggplot( )
     g <- g + ggplot2::labs( x='col', y='row')
     g <- g + ggplot2::scale_color_manual( values = getPalette(colourCount) )
     g <- g + ggplot2::geom_point(data = pgplot, ggplot2::aes(x = col, y = plot_row, color = designation), alpha = alph)
     return( g )
}


#######################################################################################
#######################################################################################
############################### SECTION: TUNABLE HELPERS ##############################
#######################################################################################
#######################################################################################

# for METRIC
TaxiNodeNodeDist <- function( node1, node2 ){
     return( abs(node1[1]-node2[1]) + abs(node1[2] - node2[2]) )
}
LinfNodeNodeDist <- function( node1, node2 ){
     return( max( abs(node1[1]-node2[1]), abs(node1[2] - node2[2]) ) )
}




# for NBHD_TENSION
TensionCalculator <- function( from_nbhd, towards_nbhd, mode = "max"){
     if( mode == "max" ){
          return( max(TENSION_MATRIX[from_nbhd,towards_nbhd], TENSION_MATRIX[towards_nbhd,from_nbhd]) )
     } else if( mode == "min" ){
          return( min(TENSION_MATRIX[from_nbhd,towards_nbhd], TENSION_MATRIX[towards_nbhd,from_nbhd]) )
     } else {
          return( TENSION_MATRIX[from_nbhd, towards_nbhd] )
     }  
}

# for TENSION_MATRIX_ENTRY
TensionEntryCalculator <- function( outward, 
                                    inward, 
                                    old_tens,
                                    min_tens,
                                    max_tens,
                                    forget_fac ){
     
     
     if( inward == 0 ){
          return( max( min_tens, forget_fac * old_tens ) )
     }
     if( inward == 1 ){
          return( min( max_tens, max( forget_fac*old_tens + .5*max_tens, old_tens ) ) )
     }
     if( inward == 2 ){
          return( min( max_tens, max( forget_fac*old_tens + .8*max_tens, old_tens) ) )
     }
     if( inward >= 3 ){
          return( max_tens )
     }
}

# for TENSION_MATRIX_ENTRY
TensionEntryCalculator_new <- function( outward, 
                                        inward, 
                                        prev_tens ){
     
     
     if( !exists('MAX_TENSION')){ 
          warning('MAX_TENSION value not found. Setting to default 2.0')
          assign('MAX_TENSION', 2, .GlobalEnv)
     }
     if( !exists('MIN_TENSION')){ 
          warning('MIN_TENSION value not found. Setting to default 0.1')
          assign('MIN_TENSION', 0.1, .GlobalEnv)
     }
     mem_offset <- MIN_TENSION
     max_min_dif <- MAX_TENSION - MIN_TENSION
     max_mem_coef <- 0.8 * max_min_dif
     
     mem_coef_fun <- function( outward, inward, max_mem_coef ){
          if( inward == 0 ){
               mem_coef <- 0
          } else {
               mem_coef <- max_mem_coef / (1 + exp(-inward))
          }
          return(mem_coef)
     }
     
     mem_persist_fun <- function( max_mem_coef, max_min_dif ){
          return( 1 - (max_mem_coef / max_min_dif) )
     }
     
     mem_coef <- mem_coef_fun( outward, inward, max_mem_coef )
     mem_persist <- mem_persist_fun( max_mem_coef, max_min_dif )
     
     tens <- mem_offset + mem_coef + (mem_persist * (prev_tens - mem_offset))
     
     return( tens )
     
}

# for DOM_MATRIX_ENTRY
DomEntryCalculator <- function( outward, 
                                inward, 
                                old_dom ){
     
     if( !exists('MAX_DOM') ){
          warning('MAX_DOM value not found. Setting to default value +Inf')
          assign('MAX_DOM', Inf, .GlobalEnv)
     }
     if( !exists('MIN_DOM') ){
          warning('MIN_DOM value not found. Setting to default value -Inf')
          assign('MIN_DOM', -Inf, .GlobalEnv)
     }
     
     dom <- min( MAX_DOM, max( MIN_DOM, outward - inward + old_dom ) )
     
     return( dom )
}

# for SHOOT
ShootDecision <- function( perp_agent, 
                           vict_agent, 
                           shoot_factor ){
     # pull out the necessary pieces
     perp_base_nbhd <- perp_agent$base_nbhd
     vict_base_nbhd <- vict_agent$base_nbhd
     perp_viol <- perp_agent$violence_fac
     # get the tension
     tens <- NBHD_TENSION( from_nbhd = perp_base_nbhd, towards_nbhd = vict_base_nbhd )
     # find the probability of shooting
     p_shoot <- 1 - exp(- shoot_factor * (perp_viol * tens)) #basically = lam * perp_viol * tens when small
     # return decision 
     return( sample( c(TRUE,FALSE), size = 1, prob = c(p_shoot,1-p_shoot) ) )
}

# for ATTACK
AttackDecision <- function( hfam, 
                            will_attack_at = 5,
                            allow_multiple_attacks = T ){
     # requires DOMINANCE MATRIX
     if( !exists('DOMINANCE_MATRIX') ) return(NULL)
     if( nrow(DOMINANCE_MATRIX) == 0 ) return(NULL)
     dom_vec <- DOMINANCE_MATRIX[ hfam, ]
     min.dom <- min( dom_vec )
     if( min.dom >= 0 ) return(NULL)
     
     dom_vec <- dom_vec[ dom_vec < 0 ]
     Tprobs <- function( x ){
          p <- min( abs(x), will_attack_at ) / will_attack_at
          return( sample( c(T,F), size = 1, prob = c(p, 1-p) ) )
     }
     to_attack <- sapply( dom_vec, FUN = Tprobs )
     if( !any( to_attack ) ) return( NULL )
     to_attack <- names( to_attack[ to_attack ] )
     if( ! allow_multiple_attacks ) to_attack <- sample( to_attack, size = 1 )
     return( to_attack )
}


############################### CONT'D: TUNABLE HELPERS ##############################
############################### SUBSECTION: PMFers ####################################


WTer_StdHostileNodeWtsByDist <- function( base_nbhd,
                                          dfun = function( x ){ x^3 / (3*log(x) + 1) } ){
     
     
     if( !PgExistChecker() ){
          warning('PLAYGROUND not found or empty. Distance weights not contributed')
          return( 0 )
     }
     
     run_anew <- function(){
          dist_from <- paste0('dist_from_',base_nbhd)
          pg <- subset(PLAYGROUND, nbhd_name != base_nbhd, select = dist_from)
          wts <- sapply( pg[,dist_from], FUN = dfun )
          names(wts) <- rownames(pg)
          return(wts)
     }
     
     if( !exists('.NODE_WTS_BY_DIST_CACHE') ){
          nwd <- list()
     } else {
          nwd <- .NODE_WTS_BY_DIST_CACHE
     }
     if( !is.null( nwd[[base_nbhd]] ) ){
          return( nwd[[base_nbhd]] )
     } else {
          nwd[[base_nbhd]] <- run_anew()
          assign( '.NODE_WTS_BY_DIST_CACHE', nwd, .GlobalEnv )
     }
     return( nwd[[base_nbhd]] )
}


WTer_AvoidEnemyNodesByDist <- function( base_nbhd,
                                        dfun = function( x ){ x^3 / (3*log(x) + 1) } ){
     
     
     if( !PgExistChecker() ){
          warning('PLAYGROUND not found. Avoid enemy weights not contributed')
          return( 0 )
     } else if( !exists('.ETHSUM') ){
          if(!CreatePlaygroundSummaries()){
               warning('PLAYGROUND summary not found. Avoid enemy weights not contributed')
               return( 0 )
          }
     }
     
     
     run_anew <- function(){
          desig <- NbhdDesignation(base_nbhd)
          dist_from <- paste0('dist_from_',desig)
          pg <- .ETHSUM[-which(names(.ETHSUM) == dist_from)]
          min_dists <- apply(pg, MARGIN=1, FUN = function(x){min(x)})
          wts <- sapply( min_dists, FUN = dfun )
          names(wts) <- rownames(pg)
          return(wts)
     }
     
     if( !exists('.NODE_WTS_BY_ENEMY_DIST_CACHE') ){
          nwd <- list()
     } else {
          nwd <- .NODE_WTS_BY_ENEMY_DIST_CACHE
     }
     if( !is.null( nwd[[base_nbhd]] ) ){
          return( nwd[[base_nbhd]] )
     } else {
          nwd[[base_nbhd]] <- run_anew()
          assign( '.NODE_WTS_BY_ENEMY_DIST_CACHE', nwd, .GlobalEnv )
     }
     return( nwd[[base_nbhd]] )
}

WTer_AvoidAuthNodes <- function( base_nbhd, 
                                 timestep, 
                                 mem_coef, 
                                 mem_persist ){
     if( !PgExistChecker() ){
          warning('PLAYGROUND not found or empty. Authority sighting weights not contributed')
          return( 0 )
     }
     if( !( exists('AUTH_SIGHTING_HIST') ) ) return( 0 )
     
     run_anew <- function( old_nbhd_list ){
          last_time <- old_nbhd_list$timestep
          old_wts <- old_nbhd_list$wts
          auth <- subset(AUTH_SIGHTING_HIST, witness_base_nbhd == base_nbhd)      
          # only allow one sighting per auth agent per base nbhd
          #auth <- dplyr::filter( .data = dplyr::group_by( .data = auth, uid, timestep ),
          #                       row_number() == sample(safe_seq_len(n()), size = 1) )
          #auth <- as.data.frame( auth )
          all_nodes <- rownames( PLAYGROUND )
          n.all_nodes <- length(all_nodes)         
          wts <- rep(0, n.all_nodes )
          names(wts) <- all_nodes
          if( is.null(last_time) | is.null(old_wts) ){
               old_wts <- wts
               last_time <- 0
               times <- unique( auth$timestep )
          } else if( timestep > last_time ) {
               times <- unique( auth$timestep[auth$timestep >= last_time] )
               "
               it is very important in 'times' to use >= in the subsetting because
               we create weights this time from the interactions last time
               "
          }
          
          for( tm in times ){
               auth_sub <- subset( auth, timestep == tm )
               time_delt <- timestep - tm
               hot_nodes <- unique( do.call( c, lapply(auth_sub$current_node, FUN = NbrIndices) ) )
               wts[hot_nodes] <- wts[hot_nodes] + (mem_coef * mem_persist^(-time_delt))
          }
          old.names <- names(old_wts)
          # don't need to multiply by mem_coef here since already has that factor from prev
          time_delt <- timestep - last_time
          wts[old.names] <- wts[old.names] + (mem_persist^(-time_delt) * old_wts)
          return(wts[wts != 0])
          
          }
     
     if( !exists('.AUTH_SIGHTING_WTS_CACHE') ){
          auth_sight <- list()
          auth_sight[[base_nbhd]] <- list()
          auth_sight[[base_nbhd]]$wts <- run_anew( auth_sight[[base_nbhd]] )
          auth_sight[[base_nbhd]]$timestep <- timestep
     } else {
          auth_sight <- .AUTH_SIGHTING_WTS_CACHE
          if( is.null(auth_sight[[base_nbhd]]) ){
               auth_sight[[base_nbhd]] <- list()
               auth_sight[[base_nbhd]]$wts <- run_anew( auth_sight[[base_nbhd]] )
               auth_sight[[base_nbhd]]$timestep <- timestep
               assign('.AUTH_SIGHTING_WTS_CACHE', auth_sight, .GlobalEnv)
          } else {
               last_time <- auth_sight[[base_nbhd]]$timestep
               if( timestep == last_time ){
                    return(auth_sight[[base_nbhd]]$wts)
               }else if( timestep > last_time ){
                    auth_sight[[base_nbhd]]$wts <- run_anew( auth_sight[[base_nbhd]] )
                    auth_sight[[base_nbhd]]$timestep <- timestep
                    assign('.AUTH_SIGHTING_WTS_CACHE', auth_sight, .GlobalEnv)
               } else {
                    warning("Error calculating auth sighting weights. No weights will be contributed.")
                    return( 0 )
               }
          }
     }
     return( auth_sight[[base_nbhd]]$wts )
}



WTer_AvoidEnemyNodes <- function( base_nbhd, 
                                 timestep, 
                                 mem_coef, 
                                 mem_persist ){
     if( !PgExistChecker() ){
          warning('PLAYGROUND not found or empty. enemyority sighting weights not contributed')
          return( 0 )
     }
     if( !( exists('ENEMY_SIGHTING_HIST') ) ) return( 0 )
     
     run_anew <- function( old_nbhd_list ){
          last_time <- old_nbhd_list$timestep
          old_wts <- old_nbhd_list$wts
          enemy <- subset(ENEMY_SIGHTING_HIST, witness_base_nbhd == base_nbhd)      
          all_nodes <- rownames( PLAYGROUND )
          n.all_nodes <- length(all_nodes)         
          wts <- rep(0, n.all_nodes )
          names(wts) <- all_nodes
          if( is.null(last_time) | is.null(old_wts) ){
               old_wts <- wts
               last_time <- 0
               times <- unique( enemy$timestep )
          } else if( timestep > last_time ) {
               times <- unique( enemy$timestep[enemy$timestep >= last_time] )
               "
               it is very important in 'times' to use >= in the subsetting because
               we create weights this time from the interactions last time
               "
          }
          for( tm in times ){
               enemy_sub <- subset( enemy, timestep == tm )
               time_delt <- timestep - tm
               hot_nodes <- unique( do.call( c, lapply(enemy_sub$observed_in_node, FUN = NbrIndices) ) )
               wts[hot_nodes] <- wts[hot_nodes] + (mem_coef * mem_persist^(-time_delt))
          }
          old.names <- names(old_wts)
          # don't need to multiply by mem_coef here since already has that factor from prev
          time_delt <- timestep - last_time
          wts[old.names] <- wts[old.names] + (mem_persist^(-time_delt) * old_wts)
          return(wts[wts != 0])
          
          }
     
     if( !exists('.ENEMY_SIGHTING_WTS_CACHE') ){
          enemy_sight <- list()
          enemy_sight[[base_nbhd]] <- list()
          enemy_sight[[base_nbhd]]$wts <- run_anew( enemy_sight[[base_nbhd]] )
          enemy_sight[[base_nbhd]]$timestep <- timestep
     } else {
          enemy_sight <- .ENEMY_SIGHTING_WTS_CACHE
          if( is.null(enemy_sight[[base_nbhd]]) ){
               enemy_sight[[base_nbhd]] <- list()
               enemy_sight[[base_nbhd]]$wts <- run_anew( enemy_sight[[base_nbhd]] )
               enemy_sight[[base_nbhd]]$timestep <- timestep
               assign('.ENEMY_SIGHTING_WTS_CACHE', enemy_sight, .GlobalEnv)
          } else {
               last_time <- enemy_sight[[base_nbhd]]$timestep
               if( timestep == last_time ){
                    return(enemy_sight[[base_nbhd]]$wts)
               }else if( timestep > last_time ){
                    enemy_sight[[base_nbhd]]$wts <- run_anew( enemy_sight[[base_nbhd]] )
                    enemy_sight[[base_nbhd]]$timestep <- timestep
                    assign('.ENEMY_SIGHTING_WTS_CACHE', enemy_sight, .GlobalEnv)
               } else {
                    warning("Error calculating enemy sighting weights. No weights will be contributed.")
                    return( 0 )
               }
          }
     }
     return( enemy_sight[[base_nbhd]]$wts )
}

WTer_AvoidHostilityNodes <- function( base_nbhd, 
                                  timestep, 
                                  mem_coef, 
                                  mem_persist,
                                  only_victim ){
     if( !PgExistChecker() ){
          warning('PLAYGROUND not found or empty. Hostility weights not contributed to hostile motion.')
          return( 0 )
     }
     if( !( exists('HOSTILITY_HIST') ) ) return( 0 )
     
     run_anew <- function( old_nbhd_list ){
          last_time <- old_nbhd_list$timestep
          old_wts <- old_nbhd_list$wts
          host_victs <- subset( HOSTILITY_HIST, (!during_attack) & (vict_base_nbhd == base_nbhd) )
          if( !only_victim ){
               host_perps <- subset( HOSTILITY_HIST, (!during_attack) & (perp_base_nbhd == base_nbhd) )
          } 
          all_nodes <- rownames( PLAYGROUND )
          n.all_nodes <- length(all_nodes)         
          wts <- rep(0, n.all_nodes )
          names(wts) <- all_nodes
          if( is.null(last_time) | is.null(old_wts) ){
               old_wts <- wts
               last_time <- 0
               times <- unique( host_victs$timestep )
               if( !only_victim ){
                    times <- unique( c(times,
                                       host_perps$timestep )
                    )
               }
          } else if( timestep > last_time ) {
               times <- unique( host_victs$timestep[host$timestep >= last_time] )
               if( !only_victim ){
                    times <- unique( c( times, 
                                        host_perps$timestep[host$timestep >= last_time] )
                    )
               }
               "
               it is very important in 'times' to use >= in the subsetting because
               we create weights this time from the interactions last time
               "
          }
          
          for( tm in times ){
               time_delt <- timestep - tm
               hot_nodes <- unique( do.call( c, lapply( host_victs[host_victs$timestep == tm, ]$vict_in_node, FUN = NbrIndices) ) )
               if( !only_victim ){
                    hot_nodes <- unique( c( hot_nodes, 
                                            do.call( c, 
                                                     lapply( host_perps[host_perps$timestep == tm, ]$perp_in_node, 
                                                             FUN = NbrIndices) ) )
                    )
               }
               wts[hot_nodes] <- wts[hot_nodes] + (mem_coef * mem_persist^(-time_delt))
          }
          old.names <- names(old_wts)
          # don't need to multiply by mem_coef here since already has that factor from prev
          time_delt <- timestep - last_time
          wts[old.names] <- wts[old.names] + (mem_persist^(-time_delt) * old_wts)
          return(wts[wts != 0])
          
          }
     
     if( !exists('.HOSTILITY_SIGHTING_WTS_CACHE') ){
          host_sight <- list()
          host_sight[[base_nbhd]] <- list()
          host_sight[[base_nbhd]]$wts <- run_anew( host_sight[[base_nbhd]] )
          host_sight[[base_nbhd]]$timestep <- timestep
     } else {
          host_sight <- .HOSTILITY_SIGHTING_WTS_CACHE
          if( is.null(host_sight[[base_nbhd]]) ){
               host_sight[[base_nbhd]] <- list()
               host_sight[[base_nbhd]]$wts <- run_anew( host_sight[[base_nbhd]] )
               host_sight[[base_nbhd]]$timestep <- timestep
               assign('.HOSTILITY_SIGHTING_WTS_CACHE', host_sight, .GlobalEnv)
          } else {
               last_time <- host_sight[[base_nbhd]]$timestep
               if( timestep == last_time ){
                    return(host_sight[[base_nbhd]]$wts)
               }else if( timestep > last_time ){
                    host_sight[[base_nbhd]]$wts <- run_anew( host_sight[[base_nbhd]] )
                    host_sight[[base_nbhd]]$timestep <- timestep
                    assign('.HOSTILITY_SIGHTING_WTS_CACHE', host_sight, .GlobalEnv)
               } else {
                    warning("Error calculating host sighting weights. No weights will be contributed.")
                    return( 0 )
               }
          }
     }
     return( host_sight[[base_nbhd]]$wts )
}


PMFer_StdHostileNodeJumper_new <- function( base_nbhd,
                                            nbhd_to, 
                                            timestep,
                                            wt_offset,
                                            adjust_by_travel_dist,
                                            adjust_by_enemy_dist,
                                            adjust_by_auth_enc,
                                            adjust_by_enemy_enc,
                                            adjust_by_hostility_enc,
                                            travel_dfun,
                                            enemy_dfun,
                                            auth_mem_coef,
                                            auth_mem_persist,
                                            enemy_mem_coef,
                                            enemy_mem_persist,
                                            hostility_mem_coef,
                                            hostility_mem_persist,
                                            only_victim ){
     if(!exists("PLAYGROUND")) stop("PLAYGROUND must exist before node jumping.")
     
     run_anew <- function(){
          if( wt_offset <= 0 ){
               warning('wt_offset must be positive for node jumping. Will jump uniformly.')
               return( NULL )
          }
          to_nodes <- rownames(subset(PLAYGROUND, nbhd_name == nbhd_to))
          n.nodes <- length(to_nodes)
          wts <- rep( 0, n.nodes )
          names(wts) <- to_nodes
          ###
          #from travel distance
          if( adjust_by_travel_dist ){
               dist_wts <- WTer_StdHostileNodeWtsByDist( base_nbhd, travel_dfun )
               names.dists <- names(dist_wts)
               match.names <- intersect( names.dists, to_nodes )
               # dist_wts is positive source since preference to small weights 
               wts[match.names] <- wts[match.names] + dist_wts[match.names]
               rm(dist_wts) # clean room from large object
          }
          ###
          #from distance to enemy designations
          if( adjust_by_enemy_dist ){
               enemy_wts <- WTer_AvoidEnemyNodesByDist( base_nbhd, enemy_dfun )
               names.enemy <- names( enemy_wts )
               match.names <- intersect( names.enemy, to_nodes )
               # enemy_wts is negative source since preference to large weights
               wts[match.names] <- wts[match.names] - enemy_wts
               rm(enemy_wts)
          }
          ###
          #from authority sightings
          if(adjust_by_auth_enc){
               auth_sight_wts <- WTer_AvoidAuthNodes( base_nbhd,
                                                      timestep, 
                                                      auth_mem_coef, 
                                                      auth_mem_persist)
               names.auth <- names(auth_sight_wts)
               match.names <- intersect(names.auth, to_nodes)
               wts[match.names] <- wts[match.names] + auth_sight_wts[match.names]
               rm(auth_sight_wts)
          }
          #from enemy sightings
          if(adjust_by_enemy_enc){
               enemy_sight_wts <- WTer_AvoidEnemyNodes( base_nbhd,
                                                      timestep, 
                                                      enemy_mem_coef, 
                                                      enemy_mem_persist)
               names.enemy <- names(enemy_sight_wts)
               match.names <- intersect(names.enemy, to_nodes)
               wts[match.names] <- wts[match.names] + enemy_sight_wts[match.names]
               rm(enemy_sight_wts)
          }
          #from hostility 
          if(adjust_by_hostility_enc){
               hostility_sight_wts <- WTer_AvoidHostilityNodes( base_nbhd,
                                                      timestep, 
                                                      hostility_mem_coef, 
                                                      hostility_mem_persist,
                                                      only_victim = only_victim)
               names.hostility <- names(hostility_sight_wts)
               match.names <- intersect(names.hostility, to_nodes)
               wts[match.names] <- wts[match.names] + hostility_sight_wts[match.names]
               rm(hostility_sight_wts)
          }
          ## shift and add offset to adjust weights
          wts <- wts + wt_offset - min(wts)
          #since these weights are for avoidance, we take the inverse
          pmf <- 1 / wts
          pmf <- pmf / sum(pmf)
          return(pmf)
     }
     
     if( !exists('.NODE_JUMPER_PMF_CACHE') ){
          njp <- list()
          njp[[base_nbhd]] <- list()
          njp[[base_nbhd]][[nbhd_to]] <- list()
          njp[[base_nbhd]][[nbhd_to]]$pmf <- run_anew()
          njp[[base_nbhd]][[nbhd_to]]$timestep <- timestep
          assign('.NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
     } else {
          njp <- .NODE_JUMPER_PMF_CACHE
          if( !is.null(njp[[base_nbhd]]) ){
               if( !is.null(njp[[base_nbhd]][[nbhd_to]]) ){
                    if( !( is.null(ts <- njp[[base_nbhd]][[nbhd_to]]$timestep) &
                           is.null(pmf <- njp[[base_nbhd]][[nbhd_to]]$pmf) ) ){
                         if( ts != timestep ){
                              njp[[base_nbhd]][[nbhd_to]]$pmf <- run_anew()
                              njp[[base_nbhd]][[nbhd_to]]$timestep <- timestep
                              assign('.NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
                         }
                    } else{
                         njp[[base_nbhd]][[nbhd_to]]$pmf <- run_anew()
                         njp[[base_nbhd]][[nbhd_to]]$timestep <- timestep
                         assign('.NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
                    }
               } else {
                    njp[[base_nbhd]][[nbhd_to]] <- list()
                    njp[[base_nbhd]][[nbhd_to]]$pmf <- run_anew()
                    njp[[base_nbhd]][[nbhd_to]]$timestep <- timestep
                    assign('.NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
               }
          } else {
               njp[[base_nbhd]] <- list()
               njp[[base_nbhd]][[nbhd_to]] <- list()
               njp[[base_nbhd]][[nbhd_to]]$pmf <- run_anew()
               njp[[base_nbhd]][[nbhd_to]]$timestep <- timestep
               assign('.NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
          }
     }
     return( njp[[base_nbhd]][[nbhd_to]]$pmf )
     
}


PMFer_StdHostileNodeJumper <- function( base_nbhd,
                                        nbhd_to, 
                                        timestep,
                                        wt_offset,
                                        travel_dfun,
                                        enemy_dfun,
                                        auth_mem_coef,
                                        auth_mem_persist ){
     if(!exists("PLAYGROUND")) stop("PLAYGROUND must exist before node jumping.")
     
     run_anew <- function(){
          if( wt_offset <= 0 ){
               warning('wt_offset must be positive for node jumping. Will jump uniformly.')
               return( NULL )
          }
          to_nodes <- rownames(subset(PLAYGROUND, nbhd_name == nbhd_to))
          n.nodes <- length(to_nodes)
          wts <- rep( 0, n.nodes )
          names(wts) <- to_nodes
          #from travel distance
          dist_wts <- WTer_StdHostileNodeWtsByDist( base_nbhd, travel_dfun )
          names.dists <- names(dist_wts)
          match.names <- intersect( names.dists, to_nodes )
          # dist_wts is positive source since preference to small weights 
          wts[match.names] <- wts[match.names] + dist_wts[match.names]
          rm(dist_wts) # clean room from large object
          ## from authority sighting
          auth_sight_wts <- WTer_AvoidAuthNodes( base_nbhd,
                                                 timestep, 
                                                 auth_mem_coef, 
                                                 auth_mem_persist)
          names.auth <- names(auth_sight_wts)
          match.names <- intersect(names.auth, to_nodes)
          # auth_sight_wts is positive source since preference to small weights 
          wts[match.names] <- wts[match.names] + auth_sight_wts[match.names]
          rm(auth_sight_wts)
          ## from distance to enemy designations
          enemy_wts <- WTer_AvoidEnemyNodesByDist( base_nbhd, enemy_dfun )
          names.enemy <- names( enemy_wts )
          match.names <- intersect( names.enemy, to_nodes )
          # enemy_wts is negative source since preference to large weights
          wts[match.names] <- wts[match.names] - enemy_wts
          rm(enemy_wts)
          ## shift and add offset to adjust weights
          wts <- wts + wt_offset - min(wts)
          #since these weights are for avoidance, we take the inverse
          pmf <- 1 / wts
          pmf <- pmf / sum(pmf)
          return(pmf)
     }
     
     if( !exists('.NODE_JUMPER_PMF_CACHE') ){
          njp <- list()
          njp[[base_nbhd]] <- list()
          njp[[base_nbhd]][[nbhd_to]] <- list()
          njp[[base_nbhd]][[nbhd_to]]$pmf <- run_anew()
          njp[[base_nbhd]][[nbhd_to]]$timestep <- timestep
          assign('.NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
     } else {
          njp <- .NODE_JUMPER_PMF_CACHE
          if( !is.null(njp[[base_nbhd]]) ){
               if( !is.null(njp[[base_nbhd]][[nbhd_to]]) ){
                    if( !( is.null(ts <- njp[[base_nbhd]][[nbhd_to]]$timestep) &
                           is.null(pmf <- njp[[base_nbhd]][[nbhd_to]]$pmf) ) ){
                         if( ts != timestep ){
                              njp[[base_nbhd]][[nbhd_to]]$pmf <- run_anew()
                              njp[[base_nbhd]][[nbhd_to]]$timestep <- timestep
                              assign('.NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
                         }
                    } else{
                         njp[[base_nbhd]][[nbhd_to]]$pmf <- run_anew()
                         njp[[base_nbhd]][[nbhd_to]]$timestep <- timestep
                         assign('.NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
                    }
               } else {
                    njp[[base_nbhd]][[nbhd_to]] <- list()
                    njp[[base_nbhd]][[nbhd_to]]$pmf <- run_anew()
                    njp[[base_nbhd]][[nbhd_to]]$timestep <- timestep
                    assign('.NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
               }
          } else {
               njp[[base_nbhd]] <- list()
               njp[[base_nbhd]][[nbhd_to]] <- list()
               njp[[base_nbhd]][[nbhd_to]]$pmf <- run_anew()
               njp[[base_nbhd]][[nbhd_to]]$timestep <- timestep
               assign('.NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
          }
     }
     return( njp[[base_nbhd]][[nbhd_to]]$pmf )
     
}


# for HOST_NBHD_JUMP_PMFer
PMFer_StdHostileNbhdJumper <- function( base_nbhd,
                                        timestep,
                                        home_wt_by_quantile, 
                                        outside_wt,
                                        ether_wt,
                                        friendly_wt,
                                        enemy_wt,
                                        dfun,
                                        ... ){
     
     if( !PgExistChecker() ) {
          warning('PLAYGROUND not found. Unable to run PMFer_StdHostileNbhdJumper.')
          return( NULL )
     }
     if( !exists('.NBHDSUM') ){
          if( !CreateNbhdSummary()) stop('.NBHDSUM not found. Unable to run PMFer_StdHostileNbhdJumper.')
     }
     
     get_home_wt <- function(){
          n.hwq <- length(home_wt_by_quantile)
          if(n.hwq < 1){
               stop("home_wt_by_quantile not found. Cannot move agents.")
          }
          else if( n.hwq == 1 ){ 
               return( home_wt_by_quantile )
          } else if( n.hwq > 1 ){
               hns <- subset(.NBHDSUM, hostile)
               n_nodes <- subset( hns, nbhd_name == base_nbhd & hostile )$n_nodes
               qs <- quantile( hns$n_nodes, probs = seq(0,1,1/n.hwq) )
               n.qs <- length( qs )
               if( n_nodes <= qs[1] ) return( home_wt_by_quantile[1] )
               if( n_nodes >= qs[n.qs]) return( home_wt_by_quantile[n.hwq] )
               for( i in 1:(n.qs-1) ){
                    if( n_nodes > qs[i] & n_nodes <= qs[i+1] ){
                         return( home_wt_by_quantile[i] )
                    }
               }
          }
          
     }
     
     run_anew <- function(){
          nbhd_wts <- list( ... )
          ifelse( is.null(match.names <- intersect( .NBHDSUM$nbhd_name, names(nbhd_wts) )),
                  nbhd_wts <- list(),
                  nbhd_wts <- nbhd_wts[match.names] )
          desig <- NbhdDesignation(base_nbhd)
          frl <- .NBHDSUM$hostile & (.NBHDSUM$designation == desig)
          enl <- .NBHDSUM$hostile & (.NBHDSUM$designation != desig)
          etl <- .NBHDSUM$designation == 'ether'
          friendly_nbhds <- unique( .NBHDSUM[frl,]$nbhd_name )
          enemy_nbhds <- unique( .NBHDSUM[enl,]$nbhd_name )
          ether_nbhds <- unique(.NBHDSUM[etl,]$nbhd_name )
          
          home_wt <- get_home_wt()
          ho_sum <- home_wt + outside_wt
          if( ho_sum <= 0 ){
               warning( 'Weights are inconsistent and unusable in PMFer_StdHostileNbhdJumper.' )
               return(NULL)
          }
          ether_wt <- max(ether_wt,0)
          friendly_wt <- max(friendly_wt,0)
          enemy_wt <- max(enemy_wt,0)
          efe_sum <- ether_wt + friendly_wt + enemy_wt
          
          dist_from_base <- paste0('dist_from_',base_nbhd)
          dists <- as.numeric(.PGSUM[,dist_from_base])
          names(dists) <- .PGSUM$nbhd_name
          nbhd_pmf <- numeric()
          if( ether_wt > 0 ){
               ether_nbhds_wts <- lapply( ether_nbhds, FUN = function(x){ 
                    if( x == base_nbhd ) return(0)
                    1.0 / dfun(dists[[x]])
               } )
               names( ether_nbhds_wts ) <- ether_nbhds
               match.names <- intersect( names(ether_nbhds_wts), names(nbhd_wts) )
               for( mn in match.names ){
                    nbhd_pmf[[mn]] <- nbhd_wts[[mn]]*ether_nbhds_wts[[mn]]
               }
               unmatch.names <- dplyr::setdiff( names(ether_nbhds_wts), names(nbhd_wts) )
               for( un in unmatch.names ){
                    nbhd_pmf[[un]] <- ether_nbhds_wts[[un]]
               }
               if( (den <- sum(nbhd_pmf[ether_nbhds])) <= 0 ){
                    # remove ether_nbhds_wts information
                    nbhd_pmf <- nbhd_pmf[ dplyr::setdiff( names(nbhd_wts), names(ether_nbhds_wts) ) ]
               } else {
                    nbhd_pmf[ ether_nbhds ] <- (outside_wt / ho_sum) * (ether_wt / efe_sum) * (nbhd_pmf[ ether_nbhds ] / den)
               }
          }
          if( friendly_wt > 0 ){
               friendly_nbhds_wts <- lapply( friendly_nbhds, FUN = function(x){ 
                    if( x == base_nbhd ) return(0)
                    1.0 / dfun(dists[[x]])
               } )
               names(friendly_nbhds_wts) <- friendly_nbhds
               match.names <- intersect( names(friendly_nbhds_wts), names(nbhd_wts) )
               for( mn in match.names ){
                    nbhd_pmf[[mn]] <- nbhd_wts[[mn]]*friendly_nbhds_wts[[mn]]
               }
               unmatch.names <- dplyr::setdiff( names(friendly_nbhds_wts), names(nbhd_wts) )
               for( un in unmatch.names ){
                    nbhd_pmf[[un]] <- friendly_nbhds_wts[[un]]
               }
               if( (den <- sum(nbhd_pmf[friendly_nbhds])) <= 0 ){
                    # remove friendly_nbhds_wts information
                    nbhd_pmf <- nbhd_pmf[ dplyr::setdiff( names(nbhd_wts), names(friendly_nbhds_wts) ) ]
               } else {
                    nbhd_pmf[ friendly_nbhds ] <- (outside_wt / ho_sum) * (friendly_wt / efe_sum) * (nbhd_pmf[ friendly_nbhds ] / den)
               }
          }
          if( enemy_wt > 0 ){
               enemy_nbhds_wts <- lapply( enemy_nbhds, FUN = function(x){ 
                    if( x == base_nbhd ) return(0)
                    1.0 / dfun(dists[[x]])
               } )
               names(enemy_nbhds_wts) <- enemy_nbhds
               match.names <- intersect( names(enemy_nbhds_wts), names(nbhd_wts) )
               for( mn in match.names ){
                    nbhd_pmf[[mn]] <- nbhd_wts[[mn]]*enemy_nbhds_wts[[mn]]
               }
               unmatch.names <- dplyr::setdiff( names(enemy_nbhds_wts), names(nbhd_wts) )
               for( un in unmatch.names ){
                    nbhd_pmf[[un]] <- enemy_nbhds_wts[[un]]
               }
               if( (den <- sum(nbhd_pmf[enemy_nbhds])) <= 0 ){
                    # remove enemy_nbhds_wts information
                    nbhd_pmf <- nbhd_pmf[ dplyr::setdiff( names(nbhd_wts), names(enemy_nbhds_wts) ) ]
               } else {
                    nbhd_pmf[ enemy_nbhds ] <- (outside_wt / ho_sum) * (enemy_wt / efe_sum) * (nbhd_pmf[ enemy_nbhds ] / den)
               }
          }
          nbhd_pmf[[base_nbhd]] <- home_wt / ho_sum
          return( nbhd_pmf )
     }
     
     
     if( ! exists('.NBHD_JUMP_PMF_CACHE') ){ 
          jump_cache <- list()
     } else {
          jump_cache <- .NBHD_JUMP_PMF_CACHE
     }
     if( is.null(jump_cache[[base_nbhd]]) ){
          pmf <- run_anew()
          jump_cache[[base_nbhd]] <- list()
          jump_cache[[base_nbhd]]$timestep <- timestep
          jump_cache[[base_nbhd]]$pmf <- pmf
          assign('.NBHD_JUMP_PMF_CACHE', jump_cache, .GlobalEnv)
     } else if( jump_cache[[base_nbhd]]$timestep != timestep ){
          pmf <- run_anew()
          jump_cache[[base_nbhd]]$timestep <- timestep
          jump_cache[[base_nbhd]]$pmf <- pmf
          assign('.NBHD_JUMP_PMF_CACHE', jump_cache, .GlobalEnv)
     } else {
          pmf <- jump_cache[[base_nbhd]]$pmf
     }
     return( pmf )
}

WTer_HostilityNodes <- function(  timestep, 
                                  mem_coef, 
                                  mem_persist ){ 
     
     if( !exists('PLAYGROUND') ) stop( 'PLAYGROUND not found. Cannot run ABM.')
     if( !exists('HOSTILITY_HIST') ) return( 0 )
     
     run_anew <- function( old_list ){
          last_time <- old_list$timestep
          old_wts <- old_list$wts
          host_hist <- subset( HOSTILITY_HIST, !during_attack )
          all_nodes <- rownames( PLAYGROUND )
          n.all_nodes <- length(all_nodes)         
          wts <- rep(0, n.all_nodes )
          names(wts) <- all_nodes
          if( is.null(last_time) | is.null(old_wts) ){
               old_wts <- wts
               last_time <- 0
               times <- unique( host_hist$timestep )
          } else if( timestep > last_time ) {
               times <- unique( host_hist$timestep[host_hist$timestep >= last_time] )
               "
               it is very important in 'times' to use >= in the subsetting because
               we create weights this time from the interactions last time
               "
          }
          
          for( tm in times ){
               hh_sub <- subset( host_hist, timestep == tm )
               time_delt <- timestep - tm
               # no unique here since there may be overlapping regions that should be awarded more
               hot_nodes <- do.call( c, lapply(hh_sub$perp_in_node, FUN = NbrIndices) )
               wts[hot_nodes] <- wts[hot_nodes] + (mem_coef * mem_persist^(-time_delt))
          }
          old.names <- names(old_wts)
          # don't need to multiply by mem_coef here since already has that factor from prev
          time_delt <- timestep - last_time
          wts[old.names] <- wts[old.names] + (mem_persist^(-time_delt) * old_wts)
          return(wts)
          
          }
     
     if( !exists('.AUTH_HOSTILITY_WTS_CACHE') ){
          hs <- list()
          hs$wts <- run_anew( hs )
          hs$timestep <- timestep
          assign('.AUTH_HOSTILITY_WTS_CACHE', hs, .GlobalEnv)
     } else {
          hs <- .AUTH_HOSTILITY_WTS_CACHE
          if( is.null(hs) ){
               hs <- list()
               hs$wts <- run_anew( hs )
               hs$timestep <- timestep
               assign('.AUTH_HOSTILITY_WTS_CACHE', hs, .GlobalEnv)
          } else {
               last_time <- hs$timestep
               if( timestep == last_time ){
                    return(hs$wts)
               }else if( timestep > last_time ){
                    hs$wts <- run_anew( hs )
                    hs$timestep <- timestep
                    assign('.AUTH_HOSTILITY_WTS_CACHE', hs, .GlobalEnv)
               } else {
                    warning("Error calculating auth sighting weights. No weights will be contributed.")
                    return( 0 )
               }
          }
     }
     
     return( hs$wts )
}



WTer_HostilityInteractionNodes <- function(  timestep, 
                                             mem_coef, 
                                             mem_persist ){ 
     
     if( !exists('PLAYGROUND') ) stop( 'PLAYGROUND not found. Cannot run ABM.')
     if( !exists('ENEMY_SIGHTING_HIST') ) return( 0 )
     
     run_anew <- function( old_list ){
          last_time <- old_list$timestep
          old_wts <- old_list$wts
          host_hist <- ENEMY_SIGHTING_HIST
          all_nodes <- rownames( PLAYGROUND )
          n.all_nodes <- length(all_nodes)         
          wts <- rep(0, n.all_nodes )
          names(wts) <- all_nodes
          if( is.null(last_time) | is.null(old_wts) ){
               old_wts <- wts
               last_time <- 0
               times <- unique( host_hist$timestep )
          } else if( timestep > last_time ) {
               times <- unique( host_hist$timestep[host_hist$timestep >= last_time] )
               "
               it is very important in 'times' to use >= in the subsetting because
               we create weights this time from the interactions last time
               "
          }
          
          for( tm in times ){
               hh_sub <- host_hist[host_hist$timestep == tm, ]
               time_delt <- timestep - tm
               # no unique outside do.call here since there may be overlapping regions that should be awarded more
               hot_nodes <- do.call( c, lapply( unique( hh_sub$witness_in_node ), FUN = NbrIndices) )
               wts[hot_nodes] <- wts[hot_nodes] + (mem_coef * mem_persist^(-time_delt))
          }
          old.names <- names(old_wts)
          # don't need to multiply by mem_coef here since already has that factor from prev
          time_delt <- timestep - last_time
          wts[old.names] <- wts[old.names] + (mem_persist^(-time_delt) * old_wts[old.names])
          return(wts)
          
          }
     
     if( !exists('.HOSTILE_INTERACTION_WTS_CACHE') ){
          hs <- list()
          hs$wts <- run_anew( hs )
          hs$timestep <- timestep
          assign('.HOSTILE_INTERACTION_WTS_CACHE', hs, .GlobalEnv)
     } else {
          hs <- .HOSTILE_INTERACTION_WTS_CACHE
          if( is.null(hs) ){
               hs <- list()
               hs$wts <- run_anew( hs )
               hs$timestep <- timestep
               assign('.HOSTILE_INTERACTION_WTS_CACHE', hs, .GlobalEnv)
          } else {
               last_time <- hs$timestep
               if( timestep == last_time ){
                    return(hs$wts)
               }else if( timestep > last_time ){
                    hs$wts <- run_anew( hs )
                    hs$timestep <- timestep
                    assign('.HOSTILE_INTERACTION_WTS_CACHE', hs, .GlobalEnv)
               } else {
                    warning("Error calculating hostile interaction weights. No weights will be contributed.")
                    return( 0 )
               }
          }
     }
     
     return( hs$wts )
     }

PMFer_AuthNodeJumper <- function(  assigned_nbhd,
                                   timestep,
                                   wt_offset,
                                   adjust_by_hostility,
                                   adjust_by_hostile_inter,
                                   hostility_mem_coef, 
                                   hostility_mem_persist,
                                   hostinter_mem_coef, 
                                   hostinter_mem_persist){
     if(!exists('PLAYGROUND')) stop("PLAYGROUND must exist before node jumping.")
     
     run_anew <- function(){
          if( wt_offset <= 0 ){
               warning('wt_offset must be positive for node jumping. Will be set to 1.')
               wt_offset <- 1
          }
          if( assigned_nbhd == 'free_'){
               ass_nodes <- rownames(PLAYGROUND)
          } else {
               if(!exists('AUTH_NBHD_NODES')){
                    warning('AUTH_NBHD_NODES list not found. Assigned authority agents will not move.')
                    return(NULL)
               }
               ass_nodes <- AUTH_NBHD_NODES[[assigned_nbhd]]
          }
          
          n.nodes <- length( ass_nodes )
          wts <- rep( 0, times = n.nodes )
          names( wts ) <- ass_nodes
          if( adjust_by_hostility ){
               host_wts <- WTer_HostilityNodes( timestep = timestep, 
                                                mem_coef = hostility_mem_coef, 
                                                mem_persist = hostility_mem_persist )
               names.host <- names(host_wts)
               match.names <- intersect( names.host, ass_nodes )
               wts[match.names] <- wts[match.names] - host_wts[match.names]
          }
          if( adjust_by_hostile_inter ){
               host_int_wts <- WTer_HostilityInteractionNodes( timestep = timestep, 
                                                               mem_coef = hostinter_mem_coef, 
                                                               mem_persist = hostinter_mem_persist )
               
               names.host <- names(host_int_wts)
               match.names <- intersect( names.host, ass_nodes )
               wts[match.names] <- wts[match.names] - host_int_wts[match.names]
          }
          
          wts <- wts + wt_offset - min(wts)
          #since these weights are for avoidance, we take the inverse
          pmf <- 1 / wts
          pmf <- pmf / sum(pmf)
          
          return( pmf )
     }
     
     if( !exists('.AUTH_NODE_JUMPER_PMF_CACHE') ){
          njp <- list()
          njp[[assigned_nbhd]] <- list()
          njp[[assigned_nbhd]]$pmf <- run_anew()
          njp[[assigned_nbhd]]$timestep <- timestep
          assign('.AUTH_NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
     } else {
          njp <- .AUTH_NODE_JUMPER_PMF_CACHE
          if( !is.null(njp[[assigned_nbhd]]) ){
               if( !is.null(ts <- njp[[assigned_nbhd]]$timestep) &
                   !is.null(pmf <- njp[[assigned_nbhd]]$pmf) ){
                    if( ts != timestep ){
                         njp[[assigned_nbhd]]$pmf <- run_anew()
                         njp[[assigned_nbhd]]$timestep <- timestep
                         assign('.AUTH_NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
                    }
               } else {
                    njp[[assigned_nbhd]]$pmf <- run_anew()
                    njp[[assigned_nbhd]]$timestep <- timestep
                    assign('.AUTH_NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
               }
          } else {
               njp[[assigned_nbhd]] <- list()
               njp[[assigned_nbhd]]$pmf <- run_anew()
               njp[[assigned_nbhd]]$timestep <- timestep
               assign('.AUTH_NODE_JUMPER_PMF_CACHE', njp, .GlobalEnv)
          }
     }
     
     return( njp[[assigned_nbhd]]$pmf )
}

AuthAgentNodeJumper <- function( auth_agents,
                                 timestep,
                                 wt_offset,
                                 adjust_by_hostility,
                                 adjust_by_hostile_inter,
                                 hostility_mem_coef, 
                                 hostility_mem_persist,
                                 hostinter_mem_coef, 
                                 hostinter_mem_persist ){
     
     #need PLAYGROUND, AUTHORITIES_DF, AUTH_NBHD_NODES
     if( is.null(auth_agents) ) return(NULL)
     if(nrow(auth_agents) == 0) return(NULL)
     assigned_nbhds <- unique( auth_agents$assigned_nbhd )
     for( an in assigned_nbhds ){
          this_an <- auth_agents$assigned_nbhd == an
          auth_sub <- auth_agents[this_an, ]
          n.ag <- nrow(auth_sub)
          if( !is.null(pmf <- PMFer_AuthNodeJumper(assigned_nbhd = an, 
                                                   timestep = timestep,
                                                   wt_offset = wt_offset,
                                                   adjust_by_hostility = adjust_by_hostility,
                                                   adjust_by_hostile_inter = adjust_by_hostile_inter,
                                                   hostility_mem_coef = hostility_mem_coef, 
                                                   hostility_mem_persist = hostility_mem_persist,
                                                   hostinter_mem_coef = hostinter_mem_coef, 
                                                   hostinter_mem_persist = hostinter_mem_persist)) ){
               auth_sub$current_node <- sample( names(pmf), size = n.ag, prob = pmf, replace = T )
          } else {
               warning(sprintf('Error creating authority pmf for nbhd %s.',an))
               auth_sub$current_node <- sample( AUTH_NBHD_NODES[[an]], size = n.ag, replace = T )
          }
          auth_agents[this_an, ]$current_node <- auth_sub$current_node
     }
     
     return( auth_agents )
}

############################### CONT'D: TUNABLE TUNABLE HELPERS #######################
############################### SUBSECTION: Movement ##################################
GoHome <- function( agents ){
     agents$current_node <- agents$base_node
     agents$current_nbhd <- agents$base_nbhd
     agents$current_nbhd_desig <- agents$designation
     return( agents )
}

NbhdDesignation <- function( nbhd_names ){
     match_list <- lapply( nbhd_names, FUN = function(x){ .NBHDSUM$designation[.NBHDSUM$nbhd_name == x] } )
     matches <- do.call(c, match_list) #Reduce( "|", match_list )
     return( matches )
}

FamDesignation <- function( fam_names ){
     match_list <- lapply( fam_names, FUN = function(x){ (.NBHDSUM$designation[.NBHDSUM$family == x])[1] } )
     matches <- do.call(c, match_list) #Reduce( "|", match_list )
     return( matches )
}

FamNbhd <- function( family ){
     return( .PGSUM$nbhd_name[.PGSUM$family == family] )
}


UpdateHostilesDF <- function( timestep ){
     # make sure HOSTILES_DF is correctly defined
     if( !HostDFExistChecker() ){
          warning('HOSTILES_DF does not exist or is empty. No movement updated.')
          return( F )
     }
     # create an alias to work on
     agents <- HOSTILES_DF
     # partition by family since attacks are family-level events
     for( fam in unique( agents$family) ){
          atk <- ATTACK( fam )
          atk_bool  <- !is.null(atk)
          # if there is an attack, update hostility records
          if( atk_bool ) UpdateAttackHostilities( fam, atk, timestep )
          in_fam <- agents$family == fam
          ag_in_fam <- agents[in_fam, ]
          # partition by nbhd; JUMPers are nbhd specific
          for( base_nbhd in unique( ag_in_fam$base_nbhd ) ){
               in_nbhd <- ag_in_fam$base_nbhd == base_nbhd
               ag_in_nbhd <- ag_in_fam[in_nbhd, ]
               # choose which nbhd to move to
               ag_in_nbhd <- HOST_NBHD_JUMPer( ag_in_nbhd, 
                                               attacking = atk_bool,
                                               timestep )
               # choose which node to move to
               ag_in_nbhd <- HOST_NODE_JUMPer( ag_in_nbhd, 
                                               attacking = atk_bool, 
                                               timestep = timestep )
               ag_in_fam[in_nbhd,] <- ag_in_nbhd
          }
          agents[in_fam, ] <- ag_in_fam
     }
     # update HOSTILES_DF
     assign('HOSTILES_DF', agents, .GlobalEnv)
     return( T )
}

"
TOM: You need to finish StdHostileNbhdJumper below and remove the excess code from UpdateHostilesDF.
Note that Jumper* class functions are the composition helpers with PMFers to be used in static functions
such as UpdateHostilesDF.
"
StdHostileNbhdJumper <- function( agents,
                                  timestep,
                                  home_wt_by_quantile, 
                                  outside_wt,
                                  ether_wt,
                                  friendly_wt,
                                  enemy_wt,
                                  dfun,
                                  adjust_for_tension,
                                  fear_fac,
                                  anger_fac,
                                  ... ){
     ### ADDED
     get_home_wt_by_quantile <- function( base_nbhd ){
          if( !adjust_for_tension ) return( home_wt_by_quantile )
          if( !(exists('TENSION_MATRIX') & exists('MIN_TENSION') & exists('MAX_TENSION')) ){
               warning( 'TENSION variables not found.' )
               return( home_wt_by_quantile )
          }
          anger <- max( TENSION_MATRIX[which(rownames(TENSION_MATRIX) == base_nbhd), ] )
          fear <- max( TENSION_MATRIX[, which(colnames(TENSION_MATRIX) == base_nbhd)] )
          fac <- 1 + fear_fac * (fear - MIN_TENSION)/(MAX_TENSION - MIN_TENSION) - 
               anger_fac * (anger - MIN_TENSION)/(MAX_TENSION - MIN_TENSION)
          return( fac * home_wt_by_quantile )
     }
     
     for( base_nbhd in unique( agents$base_nbhd ) ){
          in_nbhd <- agents$base_nbhd == base_nbhd
          
          ### ADDED
          hwbq <- get_home_wt_by_quantile( base_nbhd )
          
          pmf <- PMFer_StdHostileNbhdJumper(base_nbhd = base_nbhd,
                                            timestep = timestep,
                                            home_wt_by_quantile = home_wt_by_quantile,
                                            outside_wt = outside_wt,
                                            ether_wt = ether_wt, 
                                            friendly_wt = friendly_wt,
                                            enemy_wt = enemy_wt,
                                            dfun,
                                            ... )
          if( is.null(pmf) ){
               # trouble, go home!
               warning(sprintf("An error occured creating a nbhd jump PMF for %s at timestep %d.",base_nbhd,timestep))
               agents$current_nbhd[in_nbhd] <- base_nbhd
               agents$current_nbhd_desig[in_nbhd] <- NbhdDesignation(base_nbhd)
          } else {
               n.ag <- sum(in_nbhd)
               to_nbhds <- sample( names(pmf), size = n.ag, prob = pmf, replace = T )
               agents$current_nbhd[in_nbhd] <- to_nbhds
               agents$current_nbhd_desig[in_nbhd] <- NbhdDesignation(to_nbhds)
          }
     }
     return( agents )
}


StdHostileNodeJumper <- function( agents, 
                                  timestep,
                                  wt_offset,
                                  adjust_by_travel_dist,
                                  adjust_by_enemy_dist,
                                  adjust_by_auth_enc,
                                  adjust_by_enemy_enc,
                                  adjust_by_hostility_enc,
                                  travel_dfun,
                                  enemy_dfun,
                                  auth_mem_coef,
                                  auth_mem_persist,
                                  enemy_mem_coef,
                                  enemy_mem_persist,
                                  hostility_mem_coef,
                                  hostility_mem_persist,
                                  only_victim ){
     
     if(!exists('PLAYGROUND')) stop( "PLAYGROUND must exist before node jumping.")
     
     go_home <- agents$current_nbhd == agents$base_nbhd
     agents$current_node[go_home] <- agents$base_node[go_home]
     moved_agents <- agents[!go_home, ]
     if( nrow(moved_agents) == 0 ) return(agents)
     
     for( base_nbhd in unique(moved_agents$base_nbhd) ){
          in_nbhd <- moved_agents$base_nbhd == base_nbhd
          for( nbhd_to in unique(moved_agents$current_nbhd) ){
               to_nbhd <- moved_agents$current_nbhd == nbhd_to
               n.rem_ag <- nrow( moved_agents[in_nbhd & to_nbhd, ] )
               if( n.rem_ag <= 0 ) next
               pmf <- PMFer_StdHostileNodeJumper_new( base_nbhd = base_nbhd,
                                                      nbhd_to = nbhd_to, 
                                                      timestep = timestep,
                                                      wt_offset = wt_offset,
                                                      adjust_by_travel_dist = adjust_by_travel_dist,
                                                      adjust_by_enemy_dist = adjust_by_enemy_dist,
                                                      adjust_by_auth_enc = adjust_by_auth_enc,
                                                      adjust_by_enemy_enc = adjust_by_enemy_enc,
                                                      adjust_by_hostility_enc = adjust_by_hostility_enc,
                                                      travel_dfun = travel_dfun,
                                                      enemy_dfun = enemy_dfun,
                                                      auth_mem_coef = auth_mem_coef,
                                                      auth_mem_persist = auth_mem_persist,
                                                      enemy_mem_coef = enemy_mem_coef,
                                                      enemy_mem_persist = enemy_mem_persist,
                                                      hostility_mem_coef = hostility_mem_coef,
                                                      hostility_mem_persist = hostility_mem_persist,
                                                      only_victim = only_victim )
               if( !is.null(pmf) ){
                    moved_agents[in_nbhd & to_nbhd, ]$current_node <- sample( names(pmf), size = n.rem_ag, prob = pmf )
               } else {
                    warning( 'Error occurred creating pmf for StdHostile node jumping. Selecting uniformly.' )
                    nodes <- rownames(subset(PLAYGROUND, nbhd_name == nbhd_to))
                    if( !is.null( pmf <- PMFer_UniformByLabel(nodes) ) ){
                         moved_agents[in_nbhd & to_nbhd, ]$current_node <- sample( names(pmf), size = n.rem_ag, prob = pmf )
                    } else {
                         warning('Error in StdHostile node jumping. Sending troubled agents home.')
                         # too many errors, sending them home
                         moved_agents[in_nbhd & to_nbhd, ]$current_node <- moved_agents[in_nbhd & to_nbhd, ]$base_node
                         moved_agents[in_nbhd & to_nbhd, ]$current_nbhd <- moved_agents[in_nbhd & to_nbhd, ]$base_nbhd
                         moved_agents[in_nbhd & to_nbhd, ]$current_nbhd_desig <- moved_agents[in_nbhd & to_nbhd, ]$designation
                    }
               }
               
          }
     }
     agents[!go_home, ] <- moved_agents
     return(agents)
}


HostileMovementUpdater <- function( timestep ){
     for( fam in unique( HOSTILES_DF$family ) ){
          family <- HOSTILES_DF[HOSTILES_DF$family == fam, ]
          n.agents <- nrow( family )
          if( !is.null( attackee_fams <- ATTACK( fam ) ) ){
               # put attack code here
               ids <- PICK_ATTACKER_IDS( fam )
               which_agents <- which( family$uid %in% ids )
               attackers <- family[which_agents, ]
               family <- family[-which_agents, ]
          }
          for( ag_n in safe_seq_len(n.agents) ){
               agent <- StdHostileAgentMover( family[ag_n, ], timestep )
               HOSTILES_DF[HOSTILES_DF$uid == agent$uid, ] <- agent
               cat(sprintf('Agent %s is %s ',agent$uid,ifelse(agent$base_nbhd == agent$current_nbhd, 
                                                              'staying home.\n', 'going out.\n')))
          }
     }
     return( T )
}

########################################################################################
########################################################################################
#################################### TUNABLE VARIABLES #################################
########################################################################################
########################################################################################

.LATTICE <- RectLatticeMaker( nrows = 350, ncols = 150 ) 

HOSTILES <- list(
     folk_1 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 1, topl_c = 130, botr_r = 20, botr_c = 150 ),
          designation = 'folk',
          family = 'GD',
          agents = c( standard = (20-1+1)*(150-130+1) ) 
     ),
     folk_2 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 43, topl_c = 130, botr_r = 48, botr_c = 135 ),
          designation = 'folk',
          agents = c( standard = (48-43+1)*(135-130+1) ) 
     ),
     folk_3 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 80, topl_c = 1, botr_r = 90, botr_c = 10 ),
          designation = 'folk',
          agents = c( standard = (90-80+1)*(10-1+1) ) 
     ),
     folk_4 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 80, topl_c = 70, botr_r = 95, botr_c = 80 ),
          designation = 'folk',
          agents = c( standard = (95-80+1)*(80-70+1) ) 
     ),
     folk_5 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 100, topl_c = 30, botr_r = 105, botr_c = 40 ),
          designation = 'folk',
          agents = c( standard = (105-100+1)*(40-30+1) )
     ),
     folk_6 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 200, topl_c = 40, botr_r = 275, botr_c = 100 ),
          designation = 'folk',
          family = 'GD',
          agents = c( standard = (275-200+1)*(100-40+1) ) 
     ),
     folk_7 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 200, topl_c = 101, botr_r = 275, botr_c = 115 ),
          designation = 'folk',
          family = 'BD',
          agents = c( standard = (275-200+1)*(115-100+1) ) #one per node
     ),
     folk_8 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 325, topl_c = 40, botr_r = 350, botr_c = 60 ),
          designation = 'folk',
          family = 'GD',
          agents = c( standard = (350-325+1)*(60-40+1) ) #one per node
     ),
     folk_9 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 325, topl_c = 61, botr_r = 350, botr_c = 85 ),
          designation = 'folk',
          family = 'BD',
          agents = c( standard = (350-325+1)*(85-61+1) ) #one per node
     ),
     people_1 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 1, topl_c = 118, botr_r = 10, botr_c = 128 ),
          designation = 'people',
          agents = c( standard = (10-1+1)*(128-118+1) ) #one per node
     ),
     people_2 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 50, topl_c = 125, botr_r = 70, botr_c = 135 ),
          designation = 'people',
          family = 'LK',
          agents = c( standard = (70-50+1)*(135-125+1) ) #one per node
     ),
     people_3 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 50, topl_c = 136, botr_r = 70, botr_c = 140 ),
          designation = 'people',
          family = 'BPS',
          agents = c( standard = (70-50+1)*(140-136+1) ) #one per node
     ),
     people_4 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 50, topl_c = 141, botr_r = 70, botr_c = 145 ),
          designation = 'people',
          family = 'VL',
          agents = c( standard = (70-50+1)*(145-141+1) ) #one per node
     ),
     people_5 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 75, topl_c = 15, botr_r = 95, botr_c = 65 ),
          hostile = TRUE,
          designation = 'people',
          agents = c( standard = (95-75+1)*(65-15+1) ) #one per node
     ),
     people_6 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 200, topl_c = 5, botr_r = 215, botr_c = 35 ),
          designation = 'people',
          family = 'LK',
          agents = c( standard = (215-200+1)*(35-5+1) ) #one per node
     ),
     people_7 = list(
          nbhd_lattice = RectNbhdMaker(  topl_r = 183, topl_c = 55, botr_r = 198, botr_c = 70 ),
          designation = 'people',
          family = 'LK',
          agents = c( standard = (198-183+1)*(70-55+1) ) #one per node
     ),
     people_8 = list(
          nbhd_lattice = RectNbhdMaker(  topl_r = 183, topl_c = 71, botr_r = 198, botr_c = 105 ),
          designation = 'people',
          family = 'BPS',
          agents = c( standard = (198-183+1)*(105-71+1) ) #one per node
     ),
     people_9 = list(
          nbhd_lattice = RectNbhdMaker(  topl_r = 183, topl_c = 106, botr_r = 198, botr_c = 120 ),
          designation = 'people',
          family = 'VL',
          agents = c( standard = (198-183+1)*(120-106+1) ) #one per node
     ),
     people_10 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 277, topl_c = 40, botr_r = 323, botr_c = 90 ),
          designation = 'people',
          family = 'BPS',
          agents = c( standard = (323-277+1)*(90-40+1) ) #one per node
     ),
     people_11 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 277, topl_c = 91, botr_r = 323, botr_c = 115 ),
          designation = 'people',
          family = 'VL',
          agents = c( standard = (323-277+1)*(115-91+1) ) #one per node
     )
)

AUTHORITIES <- list( 
     n_assigned = 0,
     n_free = 500,
     assignment_radius = 10,
     min_per_nbhd = 0
)

### Any other type of nbhd
OTHER_NBHDS <- list()

### Nearby neighbor radius
NBR_RADIUS <- 2

### Use "walled neighborhoods" during agent interaction
WALLED_NBHDS <- T

### Agglomerative authority interaction
AGGL_AUTHORITY <- T

### Record all authority/hostile interactions? Or only those which prevent a hostility
ALL_AUTH_INTERACTIONS <- F

### Tension max and min
MIN_TENSION <- 0.1
MAX_TENSION <- 3

### Dominance count max and min
MAX_DOM <- Inf
MIN_DOM <- -Inf

### Run in VERBOSE mode
VERBOSE <- T

############################### CONT'D: TUNABLE VARIABLES #############################
############################### SUBSECTION: TUNABLE FUNCTIONS #########################

# METRIC must accept node1 and node2 as entries
METRIC <- function(node1, node2, ...){ TaxiNodeNodeDist(node1, node2 ) }

# NBHD_TENSION must accept from_nbhd and towards_nbhd as entries
NBHD_TENSION <- function(from_nbhd, towards_nbhd, ...){ TensionCalculator(from_nbhd, towards_nbhd, mode = 'max') }

# TENSION_MATRIX_ENTRY must accept outward, inward, and old_tens as entries
TENSION_MATRIX_ENTRY <- function( outward, inward, old_tens, ... ){ TensionEntryCalculator_new(outward,
                                                                                               inward,
                                                                                               old_tens) }

# DOM_MATRIX_ENTRY must accept outward, inward, and old_dom as entries
DOM_MATRIX_ENTRY <- function( outward, inward, old_dom, ... ){ DomEntryCalculator( outward, 
                                                                                   inward, 
                                                                                   old_dom ) }

# SHOOT must accept perp_agent and vict_agent as entries; returns T/F
SHOOT <- function( perp_agent, vict_agent, ... ){ ShootDecision(perp_agent, vict_agent, shoot_factor = 0.5) }

# ATTACK must accept attacker_fam; returns chr of family to attack, or NULL if none
ATTACK <- function( attacker_fam, ... ){  AttackDecision( attacker_fam, 
                                                          will_attack_at = 10,
                                                          allow_multiple_attacks = T ) }

ATTACK_DOM_COUNT <- function( cur_dom_cnt ){
     min_dom_cnt <- 2
     max_dom_cnt <- 6
     val <- max( min_dom_cnt, min( cur_dom_cnt + 1, max_dom_cnt ) )
     return( val )
}

HOST_NBHD_JUMPer <- function( agents, attacking, timestep, ... ){
     if( attacking ){
          # if in attack mode, how do hostile agents jump nbhds
          updated_agents <- StdHostileNbhdJumper( agents,
                                                  timestep,
                                                  home_wt_by_quantile = c(60, 70, 80, 90, 100), 
                                                  outside_wt = 2,
                                                  ether_wt = 100,
                                                  friendly_wt = 15,
                                                  enemy_wt = 0,
                                                  dfun = function(x){ dist3( dist = x, coef = 1, min_val = 1 ) },
                                                  adjust_for_tension = T,
                                                  fear_fac = .15,
                                                  anger_fac = 0.0 )
          
     } else {
          # if not in attack mode, how do hostile agents jump nbhds
          updated_agents <- StdHostileNbhdJumper( agents,
                                                  timestep,
                                                  home_wt_by_quantile = c(60, 70, 80, 90, 100), 
                                                  outside_wt = 5,
                                                  ether_wt = 100,
                                                  friendly_wt = 15,
                                                  enemy_wt = 2,
                                                  dfun = function(x){ dist3( dist = x, coef = 1, min_val = 1 ) },
                                                  adjust_for_tension = T,
                                                  fear_fac = 0.1,
                                                  anger_fac = 0.0 )
     }
     return( updated_agents )
}

HOST_NODE_JUMPer <- function( agents, attacking, timestep, ... ){ 
     # currently not using the attacking information for this stage
     updated_agents <- StdHostileNodeJumper( agents,
                                             timestep,
                                             wt_offset = 1,
                                             adjust_by_travel_dist = T,
                                             adjust_by_enemy_dist = T,
                                             adjust_by_auth_enc = F,
                                             adjust_by_enemy_enc = F,
                                             adjust_by_hostility_enc = F,
                                             travel_dfun = function(x){ 
                                                  distGraded( dist = x,
                                                              close_coef = 2,
                                                              close_cutoff = 10,
                                                              mid_cutoff = 20 )
                                             },
                                             enemy_dfun = function(x){ 
                                                  dist1( dist = x, 
                                                         coef = 1 )
                                             },
                                             auth_mem_coef = 2,
                                             auth_mem_persist = 1.5,
                                             enemy_mem_coef = 2,
                                             enemy_mem_persist = 1.5,
                                             hostility_mem_coef = 2,
                                             hostility_mem_persist = 1.5,
                                             only_victim = T  )
     return( updated_agents )
}

############################################# SETUP
.PackageHandler <- function( pkg, vers = NULL ){
     
     safe_install <- function(){
          update.packages( pkg )
          install.packages( pkg) 
     }
     
     reqs_met <- F
     if( pkg %in% installed.packages() ){
          if( !is.null(vers) ){
               if( packageVersion(pkg) < vers ){
                    cat( sprintf("Updated version of %s necessary.\n",pkg) )
                    ret <- readline( prompt = sprintf( 'Update %s now? ([yes] or no) : ', pkg) )
                    if( !( tolower(ret) == 'no' | tolower(ret) == 'n' ) ) safe_install()
               }
          }
     } else {
          ret <- readline( prompt = sprintf('Install %s now? ([y] or n): ',pkg) )
          if( !( tolower(ret) == 'no' | tolower(ret) == 'n' ) ) safe_install()
     }
     if( pkg %in% installed.packages() ) {
          if( !is.null(vers) ){
               if( packageVersion( pkg ) >= vers ) reqs_met <- T
          } else reqs_met <- T
     }
     return( reqs_met )
}

.PackagePreReqs <- function(){
     reqs_met <- F
     
     if( .PackageHandler('dplyr', vers = '0.5.0') ){
          reqs_met <- T
     } else {
          cat( paste('Package dplyr must be installed with version >= 0.5.0.\n',
                     'Please install or update this package before using this program.') )
     }
     if( !.PackageHandler('ggplot2') ) cat('Package ggplot2 is not installed. Visualization tools will fail to run.\n')
     if( !.PackageHandler('RColorBrewer') ) cat('Package RColorBrewer is not installed. Some visualization tools will fail to run.\n')
     
     return( reqs_met )
}

.TunableVariablesDefined <- function(){
     if( !exists('.LATTICE') ){
          cat(paste("ERROR: .LATTICE has not been defined.","There is no default list for .LATTICE\n", sep = " "))
          return( F )
     }
     if( !exists('HOSTILES') ){
          cat(paste("ERROR: HOSTILES list not found.","There is no default list for HOSTILES.\n", sep = " "))
          return( F )
     }
     if( !exists('AUTHORITIES') ){
          cat("- AUTHORITIES list not found. Setting to default list with 275 agents.\n")
          l <- list( n_assigned = 250, n_free = 10, assign_nbhd_by_population = T )
          assign('AUTHORITIES', l, .GlobalEnv)
     }
     if( !exists('OTHER_NBHDS') ){
          cat("-OTHER_NBHDS list not found. Setting to default empty list.\n")
          assign('OTHER_NBDHS', list(), .GlobalEnv)
     }
     if( !exists('NBR_RADIUS') ){
          cat("-NBR_RADIUS not found. Setting to default value 1.\n")
          assign('NBR_RADIUS', 1, .GlobalEnv)
     }
     if( !exists('WALLED_NBHDS') ){
          cat("-WALLED_NBHDS not found. Setting to default value TRUE.\n")
          assign('WALLED_NBHDS', T, .GlobalEnv)
     }
     if( !exists('AGGL_AUTHORITY') ){
          cat("-AGGL_AUTHORITY not found. Setting to default value TRUE.\n")
          assign('AGGL_AUTHORITY', T, .GlobalEnv)
     }
     if( !exists('VERBOSE') ){
          cat("-VERBOSE not found. Setting to default value TRUE.\n")
          assign('VERBOSE', T, .GlobalEnv)
     }
     return( T )
}

Reset <- function( ){
     rm(list = ls(pattern = '\\.[A-Z_]*_CACHE', all.names = T, envir = .GlobalEnv), envir = .GlobalEnv)
     rm(list = ls(pattern = '[A-Z_]*_HIST', envir = .GlobalEnv ), envir = .GlobalEnv )
     CreateDominanceMatrix()
     CreateTensionMatrix()
}

Setup <- function(  ){
     if(! (.PackagePreReqs()) ){
          cat("Setup failed.\n")
          return(F)
     }
     if(! (.TunableVariablesDefined()) ){
          cat("Setup failed.\n")
          return(F)
     }
     if(!CreateNbhds()){
          cat("Error creating .NBHDS list via CreateNbhds. Please check formatting of HOSTILES and OTHER_NBHDS and try again.\n")
          return(F)
     }
     if(!CreatePlayground()){ 
          cat("Error creating PLAYGROUND, .PGSUM, or .ETHSUM via CreatePlayground.\n")
          return(F)
     }
     if( !CreatePlaygroundSummaries() ){
          cat("Error creating .PGSUM and .ETHSUM via CreatePlaygroundSummaries\nSetup Failed.\n")
          return(F)
     } 
     if( !CreateNbhdSummary() ){ 
          cat("Error creating .NBHDSUM dataframe via CreateNbhdSummary\nSetup Failed.\n")
          return(F)
     } else {
          if( VERBOSE ){ 
               cat("NBHD summary:\n")
               print(.NBHDSUM )
               cat('\n')}
     }
     
     assign('HOSTILES', HostAgentCntsHelper( HOSTILES ), .GlobalEnv)
     
     if( !CreateHostilesDF() ){ 
          cat("Error creating HOSTILES_DF dataframe via CreateHostilesDF.\nSetup Failed.\n")
          return(F)
     }
     if( !CreateHostSummary() ){ 
          cat("Error creating .HOSTSUM dataframe via CreateHostSummary.\nSetup Failed.\n")
          return(F)
     } else {
          if( VERBOSE ){ 
               cat("Hostile agent summary:\n")
               print(.HOSTSUM )
               cat('\n')}
     }
     if(!CreateAuthNbhdNodes()){
          cat(paste0("Error creating AUTH_NBHD_NODES, necessasary for authority agents in model.\n",
                     "No authority agents will be added to model.\n") )
          CreateEmptyAuthoritiesDF()
     } else if(!CreateAuthoritiesDF()){ 
          cat(paste0("Error creating AUTHORITIES_DF, necessasary for authority agents in model.\n",
                     "No authority agents will be added to model.\n") )
          CreateEmptyAuthoritiesDF()
     } else if(!CreateAuthSummary()){
          cat("Error creating .HOSTNBDHSUM dataframe via CreateHostSummary.\n")
          
     } else {
          if( VERBOSE ){ 
               cat("Authority agent summary:\n")
               print(.AUTHSUM )
               cat('\n')}
     }
     if(!CreateDominanceMatrix()){
          return(F)
     }
     if( !CreateTensionMatrix() ){
          return(F)
     }
     return(T)
}

RunSample <- function( niter = 10 ){
     Reset()
     tot.s <- Sys.time()
     for( tstep in safe_seq_len(niter) ){
          t.s <- Sys.time()
          cat(sprintf('Beginning timestep %d...\n', tstep))
          cat('\t-updating hostile agents... ')
          UpdateHostilesDF( tstep )
          cat('Done.\n')
          cat('\t-updating authority agents... ')
          AUTHORITIES_DF <<- AuthAgentNodeJumper(AUTHORITIES_DF, 
                                                 tstep, 
                                                 wt_offset = 1,
                                                 adjust_by_hostility = T,
                                                 adjust_by_hostile_inter = T,
                                                 hostility_mem_coef = 100,
                                                 hostility_mem_persist = 1.25,
                                                 hostinter_mem_coef = 50,
                                                 hostinter_mem_persist = 1.25)
          cat('Done.\n')
          cat('\t-calculating interactions... ')
          UpdateInteractions(tstep)
          cat('Done.\n')
          if(exists('HOSTILITY_HIST')){
               this_time <- subset(HOSTILITY_HIST, timestep == tstep  )
               # hostilities that happen without an attack
               std_host <- dplyr::summarise(dplyr::group_by( subset(this_time, !during_attack), perp_base_nbhd, vict_base_nbhd ), 
                                            n = n())
               
               n.cnts <- nrow(std_host)
               if(n.cnts > 0){
                    cat('\t  Non-Attack Hostilities Recorded:\n')
                    for( r in safe_seq_len(nrow(std_host))){
                         inc <- std_host[r,]
                         perp <- inc$perp_base_nbhd
                         vict <- inc$vict_base_nbhd
                         cat(sprintf('\t  *%d acts with perp nbhd %s and vict nbhd %s\n', inc$n, perp, vict ))
                    }
               }
               # hostilities involved in an attack
               attack_host <- dplyr::summarise(dplyr::group_by( subset(this_time, during_attack), perp_fam, vict_fam ), 
                                               n = n())
               n.cnts <- nrow(attack_host)
               if(n.cnts > 0){
                    cat('\t Attack Hostilities Recorded:\n')
                    for( r in safe_seq_len(nrow(attack_host))){
                         inc <- attack_host[r,]
                         perp <- inc$perp_fam
                         vict <- inc$vict_fam
                         cat(sprintf('\t  *family %s attacks family %s with a domination count of %d\n', perp, vict, inc$n ))
                    }
               }
          }
          cat('\t-updating dominance and tension counts... ')
          UpdateDominanceMatrix(tstep)
          UpdateTensionMatrix(tstep)
          cat('Done.\n')
          cat('\t-updating visualization... ')
          print(ViewAgents(hsize = 1.5, asize = 1.25))
          cat('Done.\n')
          t.e <- Sys.time()
          cat( sprintf('Iteration time: %s \n\n', format(t.e - t.s)) )
     }
     tot.e <- Sys.time()
     cat(sprintf('Sample finished. Total run time: %s\n', format(tot.e-tot.s)) )
}
