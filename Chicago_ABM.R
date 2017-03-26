#################################### GENERIC HELPER FUNCTIONS
safe_seq_len <- function( n ){
     tryCatch(
          { return( seq_len( n ) ) },
          error = function(e) return( integer(length = 0L) )
     )
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

###################################### TENSION AND DOMINANCE

UpdateDominanceMatrix <- function( timestep ){
     #need DOMINANCE_MATRIX, HOSTILITY_HIST, safe_seq_len
     dm <- DOMINANCE_MATRIX
     cur_hist <- HOSTILITY_HIST[HOSTILITY_HIST$time == timestep, ]
     perp_fams <- unique( cur_hist$perp_fam )
     vict_fams <- unique( cur_hist$vict_family )
     fams <- unique( c(perp_fams, vict_fams) )
     for( i in safe_seq_len( length(fams) - 1 ) ){
          for( j in (i+1):length(fams) ){
               fi <- fams[i]
               fj <- fams[j]
               i_to_j <- sum( cur_hist$perp_family == fi & cur_hist$vict_family == fj )
               j_to_i <- sum( cur_hist$perp_family == fj & cur_hist$vict_family == fi )
               dm[fi,fj] <- DOM_MATRIX_ENTRY( outward = i_to_j, inward = j_to_i, old_dom = dm[fi,fj] )
               dm[fj,fi] <- DOM_MATRIX_ENTRY( outward = j_to_i, inward = i_to_j, old_dom = dm[fj,fi] )
          }
     }
     
     DOMINANCE_MATRIX <<- dm
     
     return(T)
}

UpdateTensionMatrix <- function( timestep ){
     #need TENSION_MATRIX, HOSTILITY_HIST, TENSION_MATRIX_ENTRY, safe_seq_len
     tm <- TENSION_MATRIX
     cur_hist <- HOSTILITY_HIST[HOSTILITY_HIST$time == timestep, ]
     perp_nbhds <- unique( cur_hist$perp_nbhd )
     vict_nbhds <- unique( cur_hist$vict_nbhd )
     nbhds <- unique( c(perp_nbhds, vict_nbhds) )
     for( i in safe_seq_len( length(nbhds) - 1 ) ){
          for( j in (i+1):length(nbhds) ){
               ni <- nbhds[i]
               nj <- nbhds[j]
               i_to_j <- sum( cur_hist$perp_family == ni & cur_hist$vict_family == nj )
               j_to_i <- sum( cur_hist$perp_family == nj & cur_hist$vict_family == ni )
               tm[ni,nj] <- TENSION_MATRIX_ENTRY( outward = i_to_j, inward = j_to_i, tm[ni,nj] )
               tm[nj,ni] <- TENSION_MATRIX_ENTRY( outward = j_to_i, inward = i_to_j, tm[nj,ni] )
          }
     }
     
     TENSION_MATRIX <- tm
     
     return(T)
}

########MOVEMENT FUNCTS


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
     
     if( exists('.NBR_INDICES_HASH') ){
          if( !is.null(.NBR_INDICES_HASH[[node_index]]) ){
               if( WALLED_NBHDS ){
                    return( .NBR_INDICES_HASH[[node_index]]$walled )
               } else {
                    return( .NBR_INDICES_HASH[[node_index]]$unwalled )
               }
          }
          nbr_indices_hash <- .NBR_INDICES_HASH
          nbr_indices_hash[[node_index]] <- list()
     } else {
          nbr_indices_hash <- list()
          nbr_indices_hash[[node_index]] <- list()
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
     nbr_indices_hash[[node_index]]$unwalled <- unwalled_nbr_idx
     pg_red <- PLAYGROUND[unwalled_nbr_idx, ]
     req_nbhd <- PLAYGROUND[node_index, ]$nbhd_name
     
     walled_nbr_idx <- rownames( pg_red[pg_red$nbhd_name == req_nbhd, ] )
     nbr_indices_hash[[node_index]]$walled <- walled_nbr_idx
     
     assign( '.NBR_INDICES_HASH', nbr_indices_hash, .GlobalEnv )
     
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
     if( ! exists('PLAYGROUND') ){
          cat('ERROR: PLAYGROUND must be created before summaries created. Try running CreatePlayground first.\n')
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
     if( !exists('PLAYGROUND') ){
          cat(paste("ERROR: Playground not found.","Has CreatePlayground been run?\n\n", sep = " "))
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

CreateHostNbhdSummary <- function( ){
     if( !exists('HOSTILES_DF') ){
          cat(paste("ERROR: HOSTILES_DF not found.",
                    "Please run CreateHostilesDF before calling this function.\n\n", sep = " "))
          return( F )
     }
     grp <- dplyr::group_by( HOSTILES_DF, designation, family, base_nbhd )
     df <- as.data.frame( dplyr::summarise( grp, n_agents = n() ) )
     assign('.HOSTNBHDSUM', df, .GlobalEnv)
     return( T )
}

CreateAuthNbhdNodes <- function( ){
     "
     Creates a list, named by neighborhoods, indicating which nodes are within
     an authority agents range of motion. This will be all nodes within the nbhd
     itself as well as all those within the ether of a certian radius specified
     within the AUTHORITIES list.
     "
     
     #make sure everything is well defined
     if( !exists('PLAYGROUND') ){
          cat(paste("ERROR: Playground not found.","Has CreatePlayground been run?\n\n", sep = " "))
          return( F )
     }
     if( !exists('AUTHORITIES') ){ 
          assign('AUTH_NBHD_NODES', list(), .GlobalEnv)
          return( T ) 
     }
     auth <- AUTHORITIES
     if( is.null(auth$assignement_radius) ) auth$assignment_radius <- 5
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
     df <- data.frame( assigned_nbhd = character(), current_node = character() )
     assign('AUTHORITIES', df, .GlobalEnv)
     return( T )
}

CreateAuthoritiesDF <- function( ){
     # Requires AUTHORITIES, AUTH_NBDH_NODES, PLAYGROUND, .PGSUM, .HOSTNBHDSUM
     
     
     if( !exists('AUTHORITIES') ){
          if( VERBOSE ){
               cat(paste("AUTHORITIES list not found.","Defaulting to no authority agents.\n\n", sep = " "))
          }
          CreateEmptyAuthoritiesDF()
          return( T )
     }
     

     auth <- AUTHORITIES
     n_free <- auth$n_free
     auth$n_free <- max( ifelse( is.null(n_free), 0, n_free ), 0 )
     if( is.null(auth$min_per_nbhd) ) auth$min_per_nbhd <- 0
     if( is.null(auth$assignement_radius) ) auth$assignment_radius <- 5
     
     
     given_nbhds <- intersect( names(AUTHORITIES), .PGSUM$nbhd_name )
     remain_nbhds <- .PGSUM$nbhd_name[ !(.PGSUM$nbhd_name %in% given_nbhds) ]
     
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
     auth$n_assigned <- max( ifelse( is.null(auth$n_assigned),0, auth$n_assigned), length.assigned )
     
     n.auth.agents.remain <- auth$n_assigned - length.assigned
     # remaining hostile nbhds 
     remain_nbhds <- .HOSTNBHDSUM$base_nbhd[ !(.HOSTNBHDSUM$base_nbhd %in% assigned_nbhds) ]
     n.nbhds.remain <- length( remain_nbhds )
     if( n.nbhds.remain == 0 ){
          auth$n_assigned <- length.assigned
          n_free <- auth$n_free + n.auth.agents.remain
     } else {
          n.auth.agents.reserved <- auth$min_per_nbhd*n.nbhds.remain
          n.auth.agents.remain <- max( n.auth.agents.remain,  n.auth.agents.reserved)
          wts <- .HOSTNBHDSUM$n_agents[ .HOSTNBHDSUM$base_nbhd %in% remain_nbhds ]
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
     
     assign('AUTHORITIES_DF', df, .GlobalEnv)
     assign('AUTHORITIES', auth, .GlobalEnv)
     
     return( T )
}

CreateAuthSummary <- function( ){
     if( !exists('AUTHORITIES_DF') ){
          cat(paste("ERROR: AUTHORITIES_DF not found.",
                    "Please run CreateAuthoritiesDF before calling this function.\n\n", sep = " "))
          return( F )
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
     tm <- SquareLabelledMatrixMaker( labs = unique( HOSTILES_DF$base_nbhd ) )
     if( is.null(tm) ){
          cat('Failed to create Tension Matrix.\n')
          return(F)
     }
     assign('TENSION_MATRIX', tm, .GlobalEnv)
     return(T)
}


############################### SUB-SECTION: INTERACTION BEHAVIOR ###################
#####################################################################################

InteractionUpdater <- function( timestep ){
     
     # retrieve all nodes neighboring an authority agent
     auth_nbrs <- NULL
     names.auth_nbrs <- NULL
     if( exists('AUTHORITIES_DF') ){
          n.auth <- nrow( AUTHORITIES_DF )
          if( n.auth > 0 ){
               #auth_nbrs <- sapply( AUTHORITIES_DF$current_node, FUN = NbrIndices )
               #auth_nbrs <- unique( do.call(c,auth_nbrs) )
               auth_nbrs <- character()
               for( i in safe_seq_len(n.auth) ){
                    auth_ag <- AUTHORITIES_DF[i, ]
                    nbrs <- NbrIndices( auth_ag$current_node )
                    for( nbr in nbrs ){
                         if( nbr %in% names.auth_nbrs ){
                              auth_nbrs[nbr] <- paste( unique( c(auth_nbrs[nbr], auth_ag$assigned_nbhd) ), collapse = "&&" )
                         } else {
                              auth_nbrs[nbr] <- auth_ag$assigned_nbhd
                              names.auth_nbrs <- names( auth_nbrs )
                         }
                    }
               }
               auth_nbrs <- auth_nbrs[ names.auth_nbrs %in% HOSTILES_DF$current_node ]
          }
     }
     # collect hostile agent info who are near an authority agent
     auth_witnesses <- HOSTILES_DF[ HOSTILES_DF$current_node %in% names.auth_nbrs, ]
     witness_id <- auth_witnesses$uid
     witness_in_node <- auth_witnesses$current_node
     witness_base_nbhd <- auth_witnesses$base_nbhd
     witness_in_nbhd <- auth_witnesses$current_nbhd
     witness_fam <- auth_witnesses$family
     witness_desig <- auth_witnesses$designation
     auth_assignment <- unname( auth_nbrs[auth_witnesses$current_node] )
     ## store these in the authority encounter log
     auth_sighting <- data.frame(  witness_id = witness_id, 
                                   witness_in_node = witness_in_node,
                                   witness_base_nbhd = witness_base_nbhd,
                                   witness_in_nbhd = witness_in_nbhd,
                                   witness_fam = witness_fam,
                                   witness_desig = witness_desig,
                                   auth_assignment = auth_assignment,
                                   stringsAsFactors = F )
     
     ## initialize the vectors that we want to remember
     ### for HOST_SIGHTING_HIST
     witness_id = character()
     witness_in_node = character()
     witness_base_nbhd = character()
     witness_in_nbhd = character()
     witness_fam = character()
     witness_desig = character()
     observed_id = character()
     observed_in_node = character()
     observed_base_nbhd = character()
     observed_fam = character()
     observed_desig = character()
     ### for HOSTILITY_HIST
     perp_id = character()
     perp_in_node = character()
     perp_nbhd = character()
     perp_fam = character()
     perp_desig = character()
     vict_id = character()
     vict_in_node = character()
     vict_nbhd = character()
     vict_fam = character()
     vict_desig = character()
     
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
               witness_id = c( witness_id, agent$uid)
               witness_in_node = c( witness_in_node, agent$current_node )
               witness_base_nbhd = c( witness_base_nbhd, agent$base_nbhd )
               witness_in_nbhd = c( witness_in_nbhd, agent$current_nbhd )
               witness_fam = c( witness_fam, agent$family )
               witness_desig = c( witness_desig, agent$designation )
               observed_id = c( observed_id, agent2$uid )
               observed_in_node = c( observed_in_node, agent2$current_node )
               observed_base_nbhd = c( observed_base_nbhd, agent2$base_nbhd )
               observed_fam = c( observed_fam, agent2$family )
               observed_desig = c( observed_desig, agent2$designation )
               if( SHOOT(perp_agent = agent, vict_agent = agent2) ){
                    # update the hostility log
                    perp_id = c( perp_id, agent$uid )
                    perp_in_node = c( perp_in_node, agent$current_node )
                    perp_nbhd = c( perp_nbhd, agent$base_nbhd )
                    perp_fam = c( perp_fam, agent$family )
                    perp_desig = c( perp_desig, agent$designation )
                    vict_id = c( vict_id, agent2$uid )
                    vict_in_node = c( vict_in_node, agent2$current_node )
                    vict_nbhd = c( vict_nbhd, agent2$base_nbhd )
                    vict_fam = c( vict_fam, agent2$family )
                    vict_desig = c( vict_desig, agent2$designation )
               }
          }
          
          who_enemies <- unmoved_agents$current_node %in% agent_nbrs
          who_enemies <- who_enemies & (unmoved_agents$designation != agent$designation)
          unmoved_enemies <- moved_agents[ who_enemies, ] 
          n.unmoved_enemies <- nrow( unmoved_enemies )
          for( mv_n2 in safe_seq_len(n.unmoved_enemies) ){
               agent2 <- unmoved_enemies[mv_n2, ]
               # update the hostility encounter log agent -> agent2
               witness_id = c( witness_id, agent$uid)
               witness_in_node = c( witness_in_node, agent$current_node )
               witness_base_nbhd = c( witness_base_nbhd, agent$base_nbhd )
               witness_in_nbhd = c( witness_in_nbhd, agent$base_nbhd )
               witness_fam = c( witness_fam, agent$family )
               witness_desig = c( witness_desig, agent$designation )
               observed_id = c( observed_id, agent2$uid )
               observed_in_node = c( observed_in_node, agent2$current_node )
               observed_base_nbhd = c( observed_base_nbhd, agent2$base_nbhd )
               observed_fam = c( observed_fam, agent2$family )
               observed_desig = c( observed_desig, agent2$designation )
               # update the hostility encounter log agent2 -> agent
               witness_id = c( witness_id, agent2$uid)
               witness_in_node = c( witness_in_node, agent2$current_node )
               witness_base_nbhd = c( witness_base_nbhd, agent2$base_nbhd )
               witness_in_nbhd = c( witness_in_nbhd, agent2$base_nbhd )
               witness_fam = c( witness_fam, agent2$family )
               witness_desig = c( witness_desig, agent2$designation )
               observed_id = c( observed_id, agent$uid )
               observed_in_node = c( observed_in_node, agent$current_node )
               observed_base_nbhd = c( observed_base_nbhd, agent$base_nbhd )
               observed_fam = c( observed_fam, agent$family )
               observed_desig = c( observed_desig, agent$designation )
               if( SHOOT(perp_agent = agent, vict_agent = agent2) ){
                    # update the hostility log agent -> agent2
                    perp_id = c( perp_id, agent$uid )
                    perp_in_node = c( perp_in_node, agent$current_node )
                    perp_nbhd = c( perp_nbhd, agent$base_nbhd )
                    perp_fam = c( perp_fam, agent$family )
                    perp_desig = c( perp_desig, agent$designation )
                    vict_id = c( vict_id, agent2$uid )
                    vict_in_node = c( vict_in_node, agent2$current_node )
                    vict_nbhd = c( vict_nbhd, agent2$base_nbhd )
                    vict_fam = c( vict_fam, agent2$family )
                    vict_desig = c( vict_desig, agent2$designation )
               }
               if( SHOOT(perp_agent = agent2, vict_agent = agent) ){
                    # update the hostility log agent2 -> agent
                    perp_id = c( perp_id, agent2$uid )
                    perp_in_node = c( perp_in_node, agent2$current_node )
                    perp_nbhd = c( perp_nbhd, agent2$base_nbhd )
                    perp_fam = c( perp_fam, agent2$family )
                    perp_desig = c( perp_desig, agent2$designation )
                    vict_id = c( vict_id, agent$uid )
                    vict_in_node = c( vict_in_node, agent$current_node )
                    vict_nbhd = c( vict_nbhd, agent$base_nbhd )
                    vict_fam = c( vict_fam, agent$family )
                    vict_desig = c( vict_desig, agent$designation )
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
                                  observed_desig = observed_desig )
     
     
     
     host_hist <- data.frame( perp_id = perp_id,
                              perp_in_node = perp_in_node,
                              perp_nbhd = perp_nbhd,
                              perp_fam = perp_fam, 
                              perp_desig = perp_desig,
                              vict_id = vict_id,
                              vict_in_node = vict_in_node,
                              vict_nbhd = vict_nbhd,
                              vict_fam = vict_fam,
                              vict_desig = vict_desig )
     
     
     
     if( nrow( host_hist ) > 0 ){
          # combine events where perps from the same nbhd, standing at the same node shoot at victs from the same nbhd
          host_hist <- dplyr::filter( .data = dplyr::group_by( .data = host_hist, perp_in_node, perp_nbhd, vict_nbhd ),
                                      row_number() == sample(safe_seq_len(n()), size = 1) )
          # combind events where multiple perps from the same nbhd shoot at the victims in the same node from the same nbhd
          host_hist <- dplyr::filter( .data = dplyr::group_by( .data = host_hist, vict_in_node, vict_nbhd, perp_nbhd ),
                                      row_number() == sample(safe_seq_len(n()), size = 1) )
          host_hist <- as.data.frame( lapply( host_hist, FUN = as.character ), stringsAsFactors = F )
     }
     ifelse( nrow(host_hist) > 0,
             host_hist$timestep <- as.integer(timestep),
             host_hist$timestep <- integer() )

     
     if( nrow( host_sighting ) > 0) {
          # combine sightings of agents from the same nbhd at the same node witnessing enemies from the same family
          host_sighting <- dplyr::filter( .data = dplyr::group_by( .data = host_sighting, witness_in_node, witness_base_nbhd, observed_fam ),
                                          row_number() == sample(safe_seq_len(n()), size = 1) )
          host_sighting <- as.data.frame( lapply( host_sighting, FUN = as.character ), stringsAsFactors = F  )
     }
     ifelse( nrow(host_sighting) > 0,
             host_sighting$timestep <- as.integer(timestep),
             host_sighting$timestep <- integer() )
     
     if( nrow( auth_sighting ) > 0 ){
          # combine sightings of agents from the same nbhd at the same node witnessing authorities
          auth_sighting <- dplyr::filter( .data = dplyr::group_by( .data = auth_sighting, witness_in_node, witness_base_nbhd ),
                                          row_number() == sample(safe_seq_len(n()), size = 1) )
          auth_sighting <- as.data.frame( lapply( auth_sighting, FUN = as.character ), stringsAsFactors = F )
     }
     ifelse( nrow(auth_sighting) > 0,
             auth_sighting$timestep <- as.integer(timestep),
             auth_sighting$timestep <- integer() )
     
     if( nrow(auth_sighting) > 0 ){
          # check if the auth_sighting prevented a hostility
          if( AGGL_AUTHORITY ){
               # record which interaction prevented hostility
               auth_to_host <- auth_sighting$witness_in_node %in% host_hist$perp_in_node
               auth_to_host <- auth_to_host | (auth_sighting$witness_in_node %in% host_hist$vict_in_node)
               auth_sighting$prevented_host <- auth_to_host
               # remove those prevented hostilities from the record (they never happened)
               host_to_auth <- host_hist$perp_in_node %in% auth_sighting$witness_in_node
               host_to_auth <- host_to_auth | (host_hist$vict_in_node %in% auth_sighting$witness_in_node)
               host_hist <- host_hist[ - which( host_to_auth ), ]
          } else {
               # record which interaction prevented hostility
               auth_to_host <- auth_sighting$witness_in_node %in% host_hist$perp_in_node
               auth_sighting$prevented_host <- auth_to_host
               # remove those prevented hostilities from the record (they never happened)
               host_to_auth <- host_hist$perp_in_node %in% auth_sighting$witness_in_node
               host_hist <- host_hist[ - which( host_to_auth ), ]
          }
     }
     if( nrow(host_hist) > 0 ){
          if( exists('HOSTILITY_HIST') ){
               host_hist <- rbind( HOSTILITY_HIST, host_hist, make.row.names = F )
          }
          assign( 'HOSTILITY_HIST', host_hist, .GlobalEnv )
     }
     if( nrow(host_sighting) > 0 ){
          if( exists('HOST_SIGHTING_HIST') ){
               host_sighting <- rbind( HOST_SIGHTING_HIST, host_sighting, make.row.names = F )
          }
          assign( 'HOST_SIGHTING_HIST', host_sighting, .GlobalEnv )
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

SavePlayground <- function( file_path ){
     saveRDS( list(.LATTICE = .LATTICE, 
                   HOSTILES = HOSTILES, 
                   PLAYGROUND = PLAYGROUND, 
                   .PGSUM = .PGSUM, 
                   .ETHSUM = .ETHSUM), 
              file = file_path )
}

OpenPlayground <- function( file_path ){
     pg <- readRDS( file_path )
     .LATTICE <<- pg$.LATTICE
     HOSTILES <<- pg$HOSTILES
     PLAYGROUND <<- pg$PLAYGROUND
     .PGSUM <<- pg$.PGSUM
     .ETHSUM <<- pg$.ETHSUM
     cat(paste0("Succesfully loaded file.\n",
               ".LATTICE, HOSTILES, PLAYGROUND, .PGSUM, .ETHSUM are now in your global environment.\n") )
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
     if(!exists('PLAYGROUND') ) stop( 'PLAYGROUND not found.' )
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
     if(!exists('PLAYGROUND') ) stop( 'PLAYGROUND not found.' )
     if( !exists('HOSTILES_DF') ) stop('HOSTILES_DF not found.')
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
     if(!exists('PLAYGROUND') ) stop( 'PLAYGROUND not found.' )
     if( !exists('AUTHORITIES_DF') ) stop('AUTHORITIES_DF not found.')
     adf <- AUTHORITIES_DF
     aplot <- .ViewAgentsPre( adf )
     
     g <- ggplot2::ggplot( aplot, ggplot2::aes(x = col, y = row) )
     g <- g + ggplot2::geom_point( color = 'firebrick3', alpha = alpha )
     return( g )
}

ViewAgents <- function( halpha = .5, aalpha = 1, only_spillover = F ){
     if(!exists('PLAYGROUND') ) stop( 'PLAYGROUND not found.' )
     if( !exists('AUTHORITIES_DF') ){
          warning('AUTHORITIES_DF not found. Only plotting hostile agents.')
          g <- ViewHostiles( alpha = halpha, only_spillover = only_spillover )
          return( g )
     }
     if( !exists('HOSTILES_DF') ) stop('HOSTILES_DF not found.')
     hdf <- HOSTILES_DF
     if( only_spillover ){
          hdf <- subset( hdf, current_nbhd != base_nbhd )
     }
     hplot <- .ViewAgentsPre( hdf )
     hplot$alpha <- 'hostile'
     hplot <- dplyr::select(hplot, col, row, family, designation, alpha )
     
     adf <- AUTHORITIES_DF
     aplot <- .ViewAgentsPre( adf )
     aplot$designation <- 'authority'
     aplot$family <- 'authority'
     aplot$alpha <- 'authority'
     aplot <- dplyr::select(aplot, col, row, family, designation, alpha)
     
     p <- rbind( hplot, aplot )
     
     
     pal <- .PalettePre( hdf )
     if( is.null(pal) ) return(NULL)
     pal <- c(pal, authority = 'firebrick3')
     
     g <- ggplot2::ggplot( p, ggplot2::aes(x = col, y = row, color = family, shape = designation, alpha = alpha ) )
     g <- g + ggplot2::scale_color_manual( values = pal )
     g <- g + ggplot2::scale_alpha_manual( values = c(hostile = halpha, authority = aalpha), guide = F )
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

# for DOM_MATRIX_ENTRY
DomEntryCalculator <- function( outward, 
                                inward, 
                                old_dom,
                                min_dom,
                                max_dom ){
     
     dom <- min( max_dom, max( min_dom, outward - inward + old_dom ) )
     
     return( dom )
}

# for SHOOT
ShootDecision <- function( perp_agent, 
                           vict_agent, 
                           shoot_factor ){
     # pull out the necessary pieces
     perp_nbhd <- perp_agent$base_nbhd
     vict_nbhd <- vict_agent$base_nbhd
     perp_viol <- perp_agent$violence_fac
     # get the tension
     tens <- NBHD_TENSION( from_nbhd = perp_nbhd, towards_nbhd = vict_nbhd )
     # find the probability of shooting
     p_shoot <- 1 - exp(- shoot_factor * (perp_viol * tens)) #basically = lam * perp_viol * tens when small
     # return decision 
     return( sample( c(TRUE,FALSE), size = 1, prob = c(p_shoot,1-p_shoot) ) )
}

# for ATTACK
AttackDecision <- function( hfam, 
                            will_attack_at = 5,
                            allow_multiple_attacks = F ){
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

# for DIST_PMFer
PMFer_ByDists <- function( dists,
                           dfun = function( x ){ x^2 },
                           offset = 0.5 ){
     dists_pmf <- 1 / (dfun(dists) + offset)
     if( (den <- sum(dists_pmf)) == 0 ) return( NULL )
     dists_pmf <- dists_pmf / den
     names( dists_pmf ) <- names( dists )
     return( dists_pmf )
}



# for AVOID_AUTH_NODES_PMFer
PMFer_HideFromAuthNodes_new <- function( base_nbhd, 
                                     nbhd_to,
                                     dist_in,
                                     nodes,
                                     timestep,
                                     mem_offset = 1,
                                     mem_coef = 1,
                                     mem_base = 2){
     
     
     diststr <- as.character(dist_in)
     
     CalculateAuthSightingNodeWts <- function( prev_wts_list ){
          # requires AUTH_SIGHTING_HIST, PLAYGROUND
          wts <- rep( mem_offset, times = length(nodes) )
          names(wts) <- nodes
          return(wts)
          if( ! exists( 'AUTH_SIGHTING_HIST' ) ) return(wts)
          
          if(is.null(prev_wts_list)){
               timedif <- 1
               prev_contr <- 0
          } else {
               if( prev_wts_list$time == timestep ) return( prev_wts_list$wts )
               timedif <- timestep - prev_wts_list$time + 1
               prev_wts <- prev_wts_list$wts - mem_offset
               prev_contr <- mem_coef * (mem_base^timedif) * prev_wts
          }
          
          auth_sighting <- subset( AUTH_SIGHTING_HIST, (witness_base_nbhd == base_nbhd) & 
                                                       (witness_in_nbhd == nbhd_to) & 
                                                       (timestep == (timestep-1) ), 
                                   select = witness_in_node )
          if( nrow(auth_sighting) == 0 ){
               wts <- wts + prev_contr[nodes]
               names(wts) <- nodes
               return( wts )
          }

          hot_nodes <- do.call( c, lapply( auth_sighting$witness_in_node, FUN = NbrIndices ) )
          hot_nodes <- intersect( hot_nodes, nodes )
          hn_sum <- dplyr::summarise(dplyr::group_by(data.frame(hn = hot_nodes), hn), n = n()) 
          
          wts <- wts + prev_contr[nodes]
          wts[hn_sum$hn] <- wts[hn_sum$hn] + ( mem_coef * hn_sum$n )
          names(wts) <- nodes
          return( wts )
     }
     
     run_anew <- function(){
          
          if( !exists('.AUTH_SIGHTING_WTS_HASH') ){
               wts <- CalculateAuthSightingNodeWts(NULL)
               asw <- list()
               asw[[base_nbhd]] <- list()
               asw[[base_nbhd]][[nbhd_to]] <- list()
               asw[[base_nbhd]][[nbhd_to]][[diststr]] <- list()
          } else {
               asw <- .AUTH_SIGHTING_WTS_HASH
               if( ! (base_nbhd %in% names(asw)) ) asw[[base_nbhd]] <- list()
               if( ! (nbhd_to %in% names(asw[[base_nbhd]])) ) asw[[base_nbhd]][[nbhd_to]] <- list()
               wts <- CalculateAuthSightingNodeWts( asw[[base_nbhd]][[nbhd_to]][[diststr]] )
          }
          asw[[base_nbhd]][[nbhd_to]][[diststr]]$time <- timestep
          asw[[base_nbhd]][[nbhd_to]][[diststr]]$wts <- wts
          assign('.AUTH_SIGHTING_WTS_HASH', asw, .GlobalEnv)
          
          if( !exists('.AUTH_SIGHTING_PMFS_HASH') ){
               ash <- list()
               ash[[base_nbhd]] <- list()
               ash[[base_nbhd]][[nbhd_to]] <- list()
               ash[[base_nbhd]][[nbhd_to]][[diststr]] <- list()
          } else {
               ash <- .AUTH_SIGHTING_PMFS_HASH
               if( ! (base_nbhd %in% names(ash)) ) ash[[base_nbhd]] <- list()
               if( ! (nbhd_to %in% names(ash[[base_nbhd]])) ) ash[[base_nbhd]][[nbhd_to]] <- list()
          }
          pmf <- 1/wts
          if( (den <- sum(pmf)) == 0 ) return(PMFer_UniformByLabel(nodes))
          pmf <- pmf / den
          ash[[base_nbhd]][[nbhd_to]][[diststr]]$time <- timestep
          ash[[base_nbhd]][[nbhd_to]][[diststr]]$pmf <- pmf
          assign('.AUTH_SIGHTING_PMFS_HASH', ash, .GlobalEnv)
          
          return( pmf )
     }
     
     if( !exists('.AUTH_SIGHTING_PMFS_HASH') ) pmf <- run_anew()
     pmf_list <- .AUTH_SIGHTING_PMFS_HASH[[base_nbhd]][[nbhd_to]][[diststr]]
     if( is.null(pmf_list) ){
          pmf <- run_anew()
     } else if( pmf_list$time != timestep ){
          pmf <- run_anew()
     } else {
          pmf <- pmf_list$pmf
     }
     
     return( pmf )
        
}


# for AVOID_AUTH_NODES_PMFer
PMFer_HideFromAuthNodes <- function( base_nbhd, 
                                     nbhd_to,
                                     dist_in,
                                     nodes,
                                     timestep,
                                     wt_offset = 1,
                                     mem_base = 2){
     
     
     
     CalculateAuthSightingNodeWts <- function( ){
          # requires AUTH_SIGHTING_HIST, PLAYGROUND
          if( ! exists( 'AUTH_SIGHTING_HIST' ) ) return(NULL)
          most_recent <- max( AUTH_SIGHTING_HIST$timestep )
          wts <- NULL
          auth <- subset(AUTH_SIGHTING_HIST, witness_base_nbhd == base_nbhd)
          for( timestep in unique(auth$timestep) ){
               cent_nodes <- auth[auth$timestep == timestep, ]$witness_in_node
               auth_sighting_nodes <- NULL
               for( node_index  in cent_nodes ){
                    auth_sighting_nodes <- c( auth_sighting_nodes, NbrIndices( node_index ) )
               }
               print(auth_sighting_nodes)
               print(nodes)
               auth_sighting_nodes <- intersect( auth_sighting_nodes, nodes )
               for( node_index in auth_sighting_nodes ){
                    if( node_index %in% names(wts) ){
                         wts[node_index] <- wts[node_index] + mem_base^{-(most_recent - timestep)}
                    } else {
                         wts <- c(wts, mem_base^{-(most_recent - timestep)})
                         names(wts)[length(wts)] <- node_index
                    }
               }
          }
          return( wts )
     }
     
     ## if we do not already have the ret. values stored, we need to recalculate them   
     
     run_anew <- function(){  
          n.nodes <- length(nodes)
          if( n.nodes == 0) return( NULL )
          if( ! exists('AUTH_SIGHTING_HIST') ) return( PMFer_UniformByLabel(nodes) )
          if( exists('AUTH_SIGHTING_HIST') & (nrow(AUTH_SIGHTING_HIST) == 0) ) return( PMFer_UniformByLabel(nodes) )
          
          auth_node_wts <- CalculateAuthSightingNodeWts()
          weighted_nodes <- intersect( nodes, names(auth_node_wts) )
          pmf <- rep( ifelse( wt_offset > 0, wt_offset, 1 ), times = n.nodes )
          names(pmf) <- nodes
          if( length(weighted_nodes) > 0 ) pmf[weighted_nodes] <- pmf[weighted_nodes] + auth_node_wts[weighted_nodes]
          if( any( pmf == 0 ) ) return( NULL )
          # move away from agents: use inverse weighting
          pmf <- 1/pmf 
          if( (den <- sum(pmf)) == 0 ) return( NULL )
          pmf <- pmf / den
          return( pmf )
     }
     
     ## use the timestep to store hashing
     timechr <- as.character(timestep)
     distchr <- as.character(dist_in)
     
     if( exists('.AUTH_SIGHTING_HASH') ){
          auth_sighting <- .AUTH_SIGHTING_HASH
          if( timechr %in% names(auth_sighting) ){
               if( base_nbhd %in% names(auth_sighting[[timechr]]) ){
                    if( nbhd_to %in% names( auth_sighting[[timechr]][[base_nbhd]]) ){
                         if( distchr %in% names( auth_sighting[[timechr]][[base_nbhd]][[distchr]] )){
                              return( auth_sighting[[timechr]][[base_nbhd]][[nbhd_to]] )
                         }
                    } else {
                         auth_sighting[[timechr]][[base_nbhd]][[nbhd_to]] <- list()
                    }
               } else {
                    auth_sighting[[timechr]][[base_nbhd]] <- list()
                    auth_sighting[[timechr]][[base_nbhd]][[nbhd_to]] <- list()
               }
          } else {
               #since we're only worried about the current timestep, can reinitialize here
               auth_sighting <- list()
               auth_sighting[[timechr]] <- list()
               auth_sighting[[timechr]][[base_nbhd]] <- list()
               auth_sighting[[timechr]][[base_nbhd]][[nbhd_to]] <- list()
          }
     } else {
          auth_sighting <- list()
          auth_sighting[[timechr]] <- list()
          auth_sighting[[timechr]][[base_nbhd]] <- list()
          auth_sighting[[timechr]][[base_nbhd]][[nbhd_to]] <- list()
     }
     
     pmf <- run_anew()
     auth_sighting[[timechr]][[base_nbhd]][[nbhd_to]][[distchr]] <- pmf
     
     ## store this so that we don't need to recalculate again
     assign('.AUTH_SIGHTING_HASH', auth_sighting, .GlobalEnv)
     
     return( pmf )
     
}


# option 1: for AVOID_ENEMY_NODES_PMFer
PMFer_HideFromEnemyNodesByDist <- function( base_nbhd,
                                            nbhd_to,
                                            dist_in,
                                            nodes,
                                            timestep,
                                            offset ){
     # currently only requires .ETHSUM
     run_anew <- function(){
          dists <- .ETHSUM[nodes, -which( names(.ETHSUM) %in% paste0('dist_from_', NbhdDesignation(base_nbhd)) ) ]
          names( dists ) <- nodes
          dists_pmf <- dists + offset
          if( (den <- sum(dists_pmf)) == 0 ) return( NULL )
          return( dists_pmf / den )
     }
     
     ## use the timestep to store hashing
     distchr <- as.character(dist_in)
     if( exists('.ENEMY_SIGHTING_HASH') ){
          enemy_hash <- .ENEMY_SIGHTING_HASH
          if( base_nbhd %in% names( enemy_hash ) ){
               if( nbhd_to %in% names( enemy_hash[[base_nbhd]]) ){
                    if( distchr %in% names( enemy_hash[[base_nbhd]][[nbhd_to]] ) ){
                         return( enemy_hash[[base_nbhd]][[nbhd_to]][[distchr]] )
                    }
               } else {
                    enemy_hash[[base_nbhd]][[nbhd_to]] <- list()
               }
          } else {
               enemy_hash[[base_nbhd]] <- list()
               enemy_hash[[base_nbhd]][[nbhd_to]] <- list()
          }
     } else {
          enemy_hash <- list()
          enemy_hash[[base_nbhd]] <- list()
          enemy_hash[[base_nbhd]][[nbhd_to]] <- list()
     }
     
     pmf <- run_anew()
     enemy_hash[[base_nbhd]][[nbhd_to]][[distchr]] <- pmf
     assign('.ENEMY_SIGHTING_HASH', enemy_hash, .GlobalEnv)
     
     return( pmf )
}


# option 2: for AVOID_ENEMY_NODES_PMFer
PMFer_HideFromHostNodesBySightings <- function( agent, 
                                                nbhd_to,
                                                dist_in,
                                                nodes,
                                                timestep,
                                                wt_offset = 1,
                                                mem_base = 2 ){
     
     
     # used in PMFer_HideFromEnemyNodesBySightings isn't based off sightings
     CalculateHostSightingNodeWts <- function( ){
          # requires HOST_SIGHTING_HIST, PLAYGROUND
          if( ! exists( 'HOST_SIGHTING_HIST' ) ) return(NULL)
          most_recent <- max( HOST_SIGHTING_HIST$timestep )
          wts <- NULL
          auth <- HOST_SIGHTING_HIST[HOST_SIGHTING_HIST$witness_base_nbhd == base_nbhd, ]
          for( timestep in unique(auth$timestep) ){
               cent_nodes <- auth[auth$timestep == timestep, ]$witness_in_node
               enemy_sighting_nodes <- NULL
               for( node_index  in cent_nodes ){
                    enemy_sighting_nodes <- c( enemy_sighting_nodes, NbrIndices( node_index ) )
               }
               enemy_sighting_nodes <- intersect( enemy_sighting_nodes, nodes )
               for( node_index in enemy_sighting_nodes ){
                    if( node_index %in% names(wts) ){
                         wts[node_index] <- wts[node_index] + mem_base^{-(most_recent - timestep)}
                    } else {
                         wts <- c(wts, mem_base^{-(most_recent - timestep)})
                         names(wts)[length(wts)] <- node_index
                    }
               }
          }
          
          return( wts )
     }
     
     ## if we do not already have the ret. values stored, we need to recalculate them   
     
     run_anew <- function(){  
          n.nodes <- length(nodes)
          if( n.nodes == 0) return( NULL )
          if( ! exists('HOST_SIGHTING_HIST') ) return( PMFer_UniformByLabel(nodes) )
          if( exists('HOST_SIGHTING_HIST') & (nrow(HOST_SIGHTING_HIST) == 0) ) return( PMFer_UniformByLabel(nodes) )
          
          host_node_wts <- CalculateHostSightingNodeWts( base_nbhd )
          weighted_nodes <- intersect( nodes, names(host_node_wts) )
          pmf <- rep( ifelse( wt_offset > 0, wt_offset, 1 ), times = n.nodes )
          names(pmf) <- nodes
          if( length(weighted_nodes) > 0 ) pmf[weighted_nodes] <- pmf[weighted_nodes] + host_node_wts[weighted_nodes]
          if( any( pmf == 0 ) ) return( NULL )
          # move away from agents: use inverse weighting
          pmf <- 1/pmf 
          if( (den <- sum(pmf)) == 0 ) return( NULL )
          pmf <- pmf / den
          return( pmf )
     }
     
     ## use the timestep to store hashing
     timechr <- as.character(timestep)
     distchr <- as.character(dist_in)
     
     if( exists('.HOST_SIGHTING_HASH') ){
          host_sighting <- .HOST_SIGHTING_HASH
          if( timechr %in% names(host_sighting) ){
               if( base_nbhd %in% names(host_sighting[[timechr]]) ){
                    if( nbhd_to %in% names( host_sighting[[timechr]][[base_nbhd]]) ){
                         if( distchr %in% names( host_sighting[[timechr]][[base_nbhd]][[distchr]] )){
                              return( host_sighting[[timechr]][[base_nbhd]][[nbhd_to]] )
                         }
                    } else {
                         host_sighting[[timechr]][[base_nbhd]][[nbhd_to]] <- list()
                    }
               } else {
                    host_sighting[[timechr]][[base_nbhd]] <- list()
                    host_sighting[[timechr]][[base_nbhd]][[nbhd_to]] <- list()
               }
          } else {
               #since we're only worried about the current timestep, can reinitialize here
               host_sighting <- list()
               host_sighting[[timechr]] <- list()
               host_sighting[[timechr]][[base_nbhd]] <- list()
               host_sighting[[timechr]][[base_nbhd]][[nbhd_to]] <- list()
          }
     } else {
          host_sighting <- list()
          host_sighting[[timechr]] <- list()
          host_sighting[[timechr]][[base_nbhd]] <- list()
          host_sighting[[timechr]][[base_nbhd]][[nbhd_to]] <- list()
     }
     
     pmf <- run_anew()
     host_sighting[[timechr]][[base_nbhd]][[nbhd_to]][[distchr]] <- pmf
     
     ## store this so that we don't need to recalculate again
     assign('.HOST_SIGHTING_HASH', host_sighting, .GlobalEnv)
     
     return( pmf )
     
}


# for HOST_NBHD_JUMP_PMFer
PMFer_StdHostMovesNbhd <- function( base_nbhd,
                                    home_wt, 
                                    ether_wt,
                                    friendly_wt,
                                    enemy_wt,
                                    ether_nbhd_wts,
                                    friendly_nbhd_wts,
                                    enemy_nbhd_wts ){
     # requires PLAYGROUND and .HOSTNBHDSUM
     make_wt_vec <- function( ){
          v <- c(home_wt = home_wt, ether_wt = ether_wt, friendly_wt = friendly_wt, enemy_wt = enemy_wt)
          if( ! ( any(v < 0 ) | (den <- sum(v) ) == 0 ) ){
               return( v/den )
          }else{
               return( NULL )
          }
     }
     make_dists <- function( desigs ){
          which_nbhds <- (PLAYGROUND$designation %in% desigs)
          # exclude the agents own nbhd
          which_nbhds <- which_nbhds & (PLAYGROUND$nbhd_name != base_nbhd)
          nbhds <- unique( PLAYGROUND$nbhd_name[which_nbhds] )
          if( length( nbhds ) == 0 ) return( NULL )
          
          dist_froms <- paste0('dist_from_',nbhds)
          red <- .PGSUM[.PGSUM$nbhd_name == base_nbhd, dist_froms ]
          dists <- as.numeric( red )
          names(dists) <- nbhds
          return( dists)
     }
     if( is.null( v <- make_wt_vec() ) ){
          stay_home <- 1
          names(stay_home) <- base_nbhd
          return( stay_home )
     }
     home_nbhd_ps <- v['home_wt']
     names(home_nbhd_ps) <- base_nbhd
     if( ether_wt > 0 ){
          if( !( any(ether_nbhd_wts < 0) | (den <- sum(ether_nbhd_wts)) == 0 ) ){
               ether_nbhd_ps <-  unname( v['ether_wt'] ) * ether_nbhd_wts / den
          } else {
               ether_nbhd_ps <- c(ether = unname( v['ether_wt']) )
          }
     } else { ether_nbhd_ps <- numeric() }
     if( friendly_wt > 0 ){
          if( !( any(friendly_nbhd_wts < 0) | (den <- sum(friendly_nbhd_wts)) == 0 ) ){
               friendly_nbhd_ps <- unname( v['friendly_wt'] ) * friendly_nbhd_wts / den
          } else {
               designation <- .HOSTNBHDSUM$designation[.HOSTNBHDSUM$base_nbhd == base_nbhd]
               dists <- make_dists( designation )
               if( is.null( dists ) ){
                    friendly_wt <- 0
                    if( is.null( v <- make_wt_vec() ) ){
                         stay_home <- 1
                         names(stay_home) <- base_nbhd
                         return( stay_home )
                    }
               } else {
                    friendly_nbhd_ps <- unname( v['friendly_wt'] ) *  DIST_PMFer( dists )
               }
          }
     } else { friendly_nbhd_ps <- numeric() }
     if( enemy_wt > 0 ){
          if( !( any(enemy_nbhd_wts < 0) | (den <- sum(enemy_nbhd_wts)) == 0 ) ){
               enemy_nbhd_ps <- unname( v['enemy_wt'] )  * enemy_nbhd_wts / den
          } else {
               enemy_desigs <- unique( .HOSTNBHDSUM$designation[.HOSTNBHDSUM$designation != designation] )
               dists <- make_dists( enemy_desigs )
               if( is.null( dists ) ){
                    enemy_wt <- 0
                    if( is.null( v <- make_wt_vec() ) ){
                         stay_home <- 1
                         names(stay_home) <- base_nbhd
                         return( stay_home )
                    }
               } else {
                    enemy_nbhd_ps <-  unname( v['enemy_wt'] ) * DIST_PMFer( dists )
               }
          }
     } else { enemy_nbhd_ps <- numeric() }
     
     pmf <- c( home_nbhd_ps, ether_nbhd_ps, friendly_nbhd_ps, enemy_nbhd_ps  )
     
     return( pmf )
}


# for NBHD_ATTACK_PMFer
PMFer_NbhdToAttack <- function( attacker_fam, 
                                attackee_fam,
                                mem_base ){
     
     # requires HOSTILITY_HIST and PLAYGROUND
     
     # nbhds belonging to the family who will be the victim of the attack
     attackee_nbhds <- unique( HOSTILES_DF$base_nbhd[ HOSTILES_DF$family == attackee_fam ] )
     print( attackee_nbhds )
     CalculateAttackNbhdWts <- function( ){
          # requires HOSTILITY_HIST, PLAYGROUND
          if( ! exists( 'HOSTILITY_HIST' ) ) return(NULL)
          if( nrow( HOSTILITY_HIST ) == 0 ) return( NULL )
          most_recent <- max( HOSTILITY_HIST$timestep )
          print( most_recent )
          wts <- NULL
          for( timestep in unique(HOSTILITY_HIST$timestep) ){
               time_host <- HOSTILITY_HIST[HOSTILITY_HIST$timestep == timestep, ]
               which_nbhds <- time_host$perp_nbhd %in% attackee_nbhds
               which_nbhds <- which_nbhds & (time_host$vict_fam == attacker_fam)
               nbhds <- time_host$perp_nbhd[ which_nbhds ]
               for( nb in nbhds ){
                    if( nb %in% names(wts) ){
                         wts[nb] <- wts[nb] + mem_base^{-(most_recent - timestep)}
                    } else {
                         wts <- c(wts, mem_base^{-(most_recent - timestep)})
                         names(wts)[length(wts)] <- nb
                    }
               }
          }
          return( wts )
     }
     
     wts <- CalculateAttackNbhdWts()
     if( is.null(wts) ) return( NULL )
     if( (den <- sum(wts)) == 0 ) return( NULL )
     pmf <- wts / den
     return( pmf )
}


# for ASSIGNED_AGENT_PMFer
PMFer_AssignedAgent <- function( agent, 
                                 timestep, 
                                 wt_offset, 
                                 mem_base ){
     # requires AUTH_NBHD_NODES
     CalculateHotNodeWts <- function( ){
          # requires HOSTILITY_HIST, PLAYGROUND
          if( ! exists( 'HOSTILITY_HIST' ) ) return(NULL)
          most_recent <- max( HOSTILITY_HIST$timestep )
          wts <- NULL
          for( timestep in unique(HOSTILITY_HIST$timestep) ){
               cent_nodes <- HOSTILITY_HIST[HOSTILITY_HIST$timestep == timestep, ]$perp_in_node
               hot_nodes <- NULL
               for( node_index  in cent_nodes ){
                    hot_nodes <- c( hot_nodes, NbrIndices( node_index ) )
               }
               for( node_index in hot_nodes ){
                    if( node_index %in% names(wts) ){
                         wts[node_index] <- wts[node_index] + mem_base^{-(most_recent - timestep)}
                    } else {
                         wts <- c(wts, mem_base^{-(most_recent - timestep)})
                         names(wts)[length(wts)] <- node_index
                    }
               }
          }
          return( wts )
     }
     run_anew <- function( ){
          hot_node_wts <- CalculateHotNodeWts()
          weighted_nodes <- intersect( AUTH_NBHD_NODES[[assigned_nbhd]], names(hot_node_wts) )
          pmf <- rep( ifelse( wt_offset > 0, wt_offset, 1 ), length(AUTH_NBHD_NODES[[assigned_nbhd]]) )
          names(pmf) <- AUTH_NBHD_NODES[[assigned_nbhd]]
          if( length(weighted_nodes) > 0 ) pmf[weighted_nodes] <- pmf[weighted_nodes] + hot_node_wts[weighted_nodes]
          if( (den <- sum(pmf)) == 0 ) return( NULL )
          pmf <- pmf / den
          return( pmf )
     }
     
     timechr <- as.character(timestep)
     assigned_nbhd <- agent$assigned_nbhd
     if( exists('.ASS_AUTH_PMF_HASH') ){
          if( !is.null( .ASS_AUTH_PMF_HASH[[timechr]] ) ){
               if( !is.null(pmf <- .ASS_AUTH_PMF_HASH[[timechr]][[assigned_nbhd]]) ){
                    return( pmf )
               } else {
                    ass_auth_hash <- .ASS_AUTH_PMF_HASH
               }
          } else {
               ass_auth_hash <- list()
               ass_auth_hash[[timechr]] <- list()
          }
     } else {
          ass_auth_hash <- list()
          ass_auth_hash[[timechr]] <- list()
     }
     
     
     pmf <- run_anew()
     ass_auth_hash[[timechr]][[assigned_nbhd]] <- pmf
     assign('.ASS_AUTH_PMF_HASH', ass_auth_hash, .GlobalEnv)
     return( pmf )
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
     match_list <- lapply( nbhd_names, FUN = function(x){ .PGSUM$designation[.PGSUM$nbhd_name == x] } )
     matches <- do.call(c, match_list) #Reduce( "|", match_list )
     return( matches )
}

StdHostileAgentsMover <- function( agents,
                                   timestep ){
     
     #go_home <- MUST_GO_HOME( agents )
     #agents[go_home, ] <- GoHome( agents[go_home, ] )
     for( base_nbhd in unique( agents$base_nbhd ) ){
          in_nbhd <- agents$base_nbhd == base_nbhd
          for( current_nbhd in unique( agents$current_nbhd[ in_nbhd ] ) ){
               in_cur_nbhd <- agents$current_nbhd == current_nbhd
               if( is.null( nbhd_pmf <- HOST_NBHD_JUMP_PMFer( base_nbhd, current_nbhd, attacking = F, timestep = timestep ) ) ){
                    warning( paste("Error creating PMF for standard hostile agent move at time:", timestep, sep = " ") )
                    agents$current_nbhd[ in_nbhd & in_cur_nbhd ] <- base_nbhd
                    agents$current_nbhd_desig[ in_nbhd & in_cur_nbhd ] <- NbhdDesignation( base_nbhd )
               } else {
                    n.ag <- sum(in_nbhd & in_cur_nbhd)
                    nbhds_to <- sample( names(nbhd_pmf), size = n.ag, replace = T, prob = nbhd_pmf )
                    agents$current_nbhd[ in_nbhd & in_cur_nbhd ] <- nbhds_to
                    #return( agents[ in_nbhd & in_cur_nbhd, ])
                    agents$current_nbhd_desig[ in_nbhd & in_cur_nbhd ] <- NbhdDesignation( nbhds_to )
               }
          }
          "
          At this point, what we have left for this iteration is to find the node each agent jumps to
          based on the updated current_nbhd (and current_nbhd_desig) values. These agents are passed 
          to the global function HOST_NODE_JUMPer, which returns a list of nodes.
          "
          agents$current_node[ in_nbhd ] <- HOST_NODE_JUMPer( agents[ in_nbhd, ], 
                                                                           attacking = F, 
                                                                           timestep = timestep )
          
     }
     
     return( agents )
     
}

StdHostileNodeJumper <- function( agents, 
                                  timestep,
                                  auth_offset,
                                  auth_mem_base, 
                                  enemy_dist_offset,
                                  auth_wt,
                                  enemy_wt ){
     go_home <- agents$current_nbhd == agents$base_nbhd
     agents$current_node[go_home] <- agents$base_node[go_home]
     moved_agents <- agents[!go_home, ]
     if( nrow(moved_agents) == 0 ) return(agents)
     for( base_nbhd in unique(moved_agents$base_nbhd) ){
          in_nbhd <- moved_agents$base_nbhd == base_nbhd
          for( current_nbhd in unique( moved_agents$current_nbhd[ in_nbhd ] ) ){
               in_cur_nbhd <- moved_agents$current_nbhd == current_nbhd
               n.ag <- nrow( moved_agents[ in_nbhd & in_cur_nbhd, ] )
               pg <- subset(PLAYGROUND, nbhd_name == current_nbhd)
               dist_from <- paste0('dist_from_', base_nbhd)
               dists <- pg[ , dist_from]
               if( is.null( dists_pmf <- DIST_PMFer(dists) ) ){
                    possible_nodes <- rownames( pg )
                    moved_agents$current_node[ in_nbhd & in_cur_nbhd ] <- sample( possible_nodes, size = n.ag, replace = T)
               } else {
                    dists_in <- sample( dists, size = n.ag, replace = T, prob = dists_pmf )
                    for( dist_in in unique(dists_in) ){
                         at_dist <- dists_in == dist_in
                         n.at_dist <- sum(at_dist)
                         pg_matches_dist <- pg[,dist_from] == dist_in
                         nodes <- rownames( pg[pg_matches_dist, ] )
                         avoid_enemy_pmf <- NULL
                         avoid_auth_pmf <- NULL
                         avoid_enemy_pmf <- PMFer_HideFromEnemyNodesByDist(base_nbhd, 
                                                                           current_nbhd, 
                                                                           dist_in, 
                                                                           nodes, 
                                                                           timestep, 
                                                                           offset = enemy_dist_offset)
                         if( is.null( avoid_enemy_pmf ) ){
                              avoid_enemy_pmf <- PMFer_UniformByLabel( nodes )
                         }
                         avoid_auth_pmf <- PMFer_HideFromAuthNodes_new( base_nbhd, 
                                                                      current_nbhd,
                                                                      dist_in,
                                                                      nodes,
                                                                      timestep,
                                                                      mem_offset = auth_offset,
                                                                      mem_base = auth_mem_base )
                         if( is.null( avoid_auth_pmf ) ){
                              avoid_auth_pmf <- PMFer_UniformByLabel( nodes )
                         }
                         
                         ## create a weighted average of the avoid_auth/enemy_pmfs
                         pmf <- PMFer_CombineWeightedPMFs( wts = c(auth_wt, enemy_wt), avoid_auth_pmf, avoid_enemy_pmf )
                         moved_agents$current_node[ in_nbhd & in_cur_nbhd & at_dist ] <- sample( nodes, size = n.at_dist, replace = T, prob = pmf )
                    }
               }
          }
     }
     agents[!go_home, ] <- moved_agents
     return( agents )
}

HOST_NODE_JUMPer <- function( agents, attacking, timestep, ... ){ 
     if( !attacking ){
          updated_agents <- StdHostileNodeJumper( agents, 
                                                    timestep,
                                                    auth_offset = 1,
                                                    auth_mem_base = 2, 
                                                    enemy_dist_offset = 5,
                                                    auth_wt = 1,
                                                    enemy_wt = 1 )
          return( updated_agents$current_node )
     } else {
          AttackingHostileNodeJumper( agents )
     }
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
          agents = c( standard = (20-1+1)*(150-130+1) ) #one per node
     ),
     folk_2 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 43, topl_c = 130, botr_r = 48, botr_c = 135 ),
          designation = 'folk',
          agents = c( standard = (48-43+1)*(135-130+1) ) #one per node
     ),
     folk_3 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 80, topl_c = 1, botr_r = 90, botr_c = 10 ),
          designation = 'folk',
          agents = c( standard = (90-80+1)*(10-1+1) ) #one per node
     ),
     folk_4 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 80, topl_c = 70, botr_r = 95, botr_c = 80 ),
          designation = 'folk',
          agents = c( standard = (95-80+1)*(80-70+1) ) #one per node
     ),
     folk_5 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 100, topl_c = 30, botr_r = 105, botr_c = 40 ),
          designation = 'folk',
          agents = c( standard = (105-100+1)*(40-30+1) ) #one per node
     ),
     folk_6 = list(
          nbhd_lattice = RectNbhdMaker( topl_r = 200, topl_c = 40, botr_r = 275, botr_c = 100 ),
          designation = 'folk',
          family = 'GD',
          agents = c( standard = (275-200+1)*(100-40+1) ) #one per node
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
     n_assigned = 120,
     n_free = 10,
     assignment_radius = 5,
     min_per_nbhd = 3
)

### Any other type of nbhd
OTHER_NBHDS <- list()

### Nearby neighbor radius
NBR_RADIUS <- 1

### Use "walled neighborhoods" during agent interaction
WALLED_NBHDS <- T

### Agglomerative authority interaction
AGGL_AUTHORITY <- T

### Run in VERBOSE mode
VERBOSE <- T

############################### CONT'D: TUNABLE VARIABLES #############################
############################### SUBSECTION: TUNABLE FUNCTIONS #########################

# METRIC must accept node1 and node2 as entries
METRIC <- function(node1, node2, ...){ TaxiNodeNodeDist(node1, node2 ) }

# NBHD_TENSION must accept from_nbhd and towards_nbhd as entries
NBHD_TENSION <- function(from_nbhd, towards_nbhd, ...){ TensionCalculator(from_nbhd, towards_nbhd, mode = 'max') }

# TENSION_MATRIX_ENTRY must accept outward, inward, and old_tens as entries
TENSION_MATRIX_ENTRY <- function( outward, inward, old_tens, ... ){ TensionEntryCalculator(outward,
                                                                                           inward,
                                                                                           old_tens,
                                                                                           min_tens = 0.1,
                                                                                           max_tens = 1,
                                                                                           forget_fac = 0.95) }
# DOM_MATRIX_ENTRY must accept outward, inward, and old_dom as entries
DOM_MATRIX_ENTRY <- function( outward, inward, old_dom, ... ){ DomEntryCalculator( outward, 
                                                                                   inward, 
                                                                                   old_dom, 
                                                                                   min_dom = -Inf,
                                                                                   max_dom = Inf) }

# SHOOT must accept perp_agent and vict_agent as entries; returns T/F
SHOOT <- function( perp_agent, vict_agent, ... ){ ShootDecision(perp_agent, vict_agent, shoot_factor = 0.1) }

# ATTACK must accept attacker_fam; returns chr of family to attack, or NULL if none
ATTACK <- function( attacker_fam, ... ){ AttackDecision( attacker_fam, 
                                                           will_attack_at = 5,
                                                           allow_multiple_attacks = F ) }

# N_ATTACKERS must accept attacker_fam and attackee_fams
PICK_ATTACKER_IDS <- function( attacker_fam, attackee_fams, ... ){ NULL }

### PMFers and movement

# MUST_GO_HOME
# MUST_GO_HOME <- function( agents, ... ){ return( agents$current_nbhd != agents$base_nbhd ) }


# ASSIGNED_AGENT_PMFer must accept agent and timestamp
ASSIGNED_AGENT_PMFer <- function( agent, timestep, ... ){ PMFer_AssignedAgent( agent, 
                                                                               timestep, 
                                                                               wt_offset = 1, 
                                                                               mem_base = 2 ) }
# NBHD_ATTACK_PMFer must accept attacker_fam and attackee_fam
NBHD_ATTACK_PMFer <- function( attacker_fam, attackee_fam, ... ){ PMFer_NbhdToAttack( attacker_fam, 
                                                                                      attackee_fam,
                                                                                      mem_base = 1.5 ) }

# DIST_PMFer must accept dists
DIST_PMFer <- function( dists, ... ){ PMFer_ByDists(dists, 
                                                    dfun = function(x){ x^2},
                                                    offset = 0.5 ) }

# HOST_NBHD_JUMP_PMFer must accept agent, attacking (T/F whether an attack is taking place) and timestep
HOST_NBHD_JUMP_PMFer <- function( base_nbhd, current_nbhd, attacking, timestep, ... ){
     if( attacking ){
          home_wt <- .95
          ether_wt <- .05
          friendly_wt <- 0 
          enemy_wt <- 0
     } else {
          home_wt <- .85
          ether_wt <- .15
          friendly_wt <- 0 
          enemy_wt <- 0
     }
     PMFer_StdHostMovesNbhd( base_nbhd,
                             home_wt, 
                             ether_wt,
                             friendly_wt,
                             enemy_wt,
                             ether_nbhd_wts = NULL,
                             friendly_nbhd_wts = NULL,
                             enemy_nbhd_wts = NULL )
}


# HIDE_FROM_ENEMY_NODES_PMFer
AVOID_ENEMY_NODES_PMFer <- function( agent, nbhd_to, dist_in, nodes, attacking, timestep, ...){ 
     PMFer_HideFromEnemyNodesByDist( agent,
                                     nbhd_to,
                                     dist_in,
                                     nodes,
                                     timestep,
                                     offset = 5 )
}

# AVOID_AUTH_NODES_PMFer
AVOID_AUTH_NODES_PMFer <- function( agent, nbhd_to, dist_in, nodes, attacking, timestep, ...){ 
     PMFer_HideFromAuthNodes( agent, 
                              nbhd_to,
                              dist_in,
                              nodes,
                              timestep,
                              wt_offset = 1,
                              mem_base = 2)
}

# HOST_NODE_PMFer
HOST_NODE_JUMP_PMFer <- function( agent, nbhd_to, attacking, timestep, ...){
     ## go home if staying in the same nbhd
     if( nbhd_to == agent$base_nbhd ){
          ret <- 1
          names(ret) <- agent$base_node
          return( ret ) 
     }
     ## distance to travel into nbhd_to (tunable via DIST_PMFer)
     dists <- unique( PLAYGROUND[PLAYGROUND$nbhd_name == nbhd_to, paste0('dist_from_', agent$base_nbhd)] )
     dists_pmf <- DIST_PMFer( dists )
     ## go home if there is an error
     if( is.null( dists_pmf ) ){
          ret <- 1
          names(ret) <- agent$base_node
          return( ret ) 
     }
     dist_in <- sample( dists, size = 1, prob = dists_pmf )
     
     pg <- PLAYGROUND[ PLAYGROUND$nbhd_name == nbhd_to, ]
     pg <- pg[ pg[[paste0('dist_from_',agent$base_nbhd)]] == dist_in, ]
     nodes <- rownames( pg )
     ## get weighted pmfs based on authority interaction and enemy locations/interactions
     avoid_enemy_pmf <- AVOID_ENEMY_NODES_PMFer(agent, nbhd_to, dist_in, nodes, attacking, timestep, ...)
     if( is.null( avoid_enemy_pmf ) ){
          avoid_enemy_pmf <- PMFer_UniformByLabel( nodes )
     }
     avoid_auth_pmf <- AVOID_AUTH_NODES_PMFer(agent, nbhd_to, dist_in, nodes, attacking, timestep, ...)
     if( is.null( avoid_auth_pmf ) ){
          avoid_auth_pmf <- PMFer_UniformByLabel( nodes )
     }
     
     ## create a weighted average of the avoid_auth/enemy_pmfs
     ## tunable here via weights assigned to each pmf
     if( !attacking ){
          # equal weighting to avoid auth and enemy
          auth_wt <- 1
          enemy_wt <- 1
     } else {
          # more likely to avoid enemy 
          auth_wt <- 1
          enemy_wt <- 2
     }
     PMFer_CombineWeightedPMFs( wts = c(enemy_wt, auth_wt), avoid_enemy_pmf, avoid_auth_pmf )
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
     if( !CreateHostilesDF() ){ 
          cat("Error creating HOSTILES_DF dataframe via CreateHostilesDF.\nSetup Failed.\n")
          return(F)
     }
     if( !CreateHostNbhdSummary() ){ 
          cat("Error creating .HOSTNBDHSUM dataframe via CreateHostNbhdSummary.\nSetup Failed.\n")
          return(F)
     } else {
          if( VERBOSE ){ 
               cat("Hostile agent summary:\n")
               print(.HOSTNBHDSUM )
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
          cat("Error creating .HOSTNBDHSUM dataframe via CreateHostNbhdSummary.\n")
         
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
