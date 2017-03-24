#################################### GENERIC HELPER FUNCTIONS
seq_maker <- function( n ){
     if( class(n) == "numeric" ){ n <- as.integer(n) }
     if( class(n) != "integer" ) return(NULL)
     return( seq_len( max(n,0) ) )
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
     return( pmf/den )
}

###################################### TENSION AND DOMINANCE

UpdateDominanceMatrix <- function( timestep ){
     #need DOMINANCE_MATRIX, HOSTILITY_HIST, seq_maker
     dm <- DOMINANCE_MATRIX
     cur_hist <- HOSTILITY_HIST[HOSTILITY_HIST$time == timestep, ]
     perp_fams <- unique( cur_hist$perp_fam )
     vict_fams <- unique( cur_hist$vict_family )
     fams <- unique( c(perp_fams, vict_fams) )
     for( i in seq_maker( length(fams) - 1 ) ){
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
     #need TENSION_MATRIX, HOSTILITY_HIST, TENSION_MATRIX_ENTRY, seq_maker
     tm <- TENSION_MATRIX
     cur_hist <- HOSTILITY_HIST[HOSTILITY_HIST$time == timestep, ]
     perp_nbhds <- unique( cur_hist$perp_nbhd )
     vict_nbhds <- unique( cur_hist$vict_nbhd )
     nbhds <- unique( c(perp_nbhds, vict_nbhds) )
     for( i in seq_maker( length(nbhds) - 1 ) ){
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
     playground <- merge( playground, LATTICE[,c('row','col')], by = c('row','col'), all.x = T, all.y = T)
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
     for( an in seq_maker(n.non_ether_nbhds) ){
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
     nb_grp <- dplyr::group_by(PLAYGROUND, nbhd_name )
     pgsum <-   dplyr::summarise_at(nb_grp, .cols = dplyr::vars(dplyr::contains("dist_from_")), .funs = min )
     pgsum$hostile <- as.logical(
          ( dplyr::summarise(nb_grp, h = max(hostile)) )$h
     )
     assign('.PGSUM', pgsum, .GlobalEnv )
     
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
     assign('.HOSTNBHDSUM', dplyr::summarise( grp, n_agents = n() ), .GlobalEnv)
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
     auth$n_free <- max( ifelse( is.null(auth$n_free), 0, auth$n_free ), 0 )
     if( is.null(auth$min_per_nbhd) ) auth$min_per_nbhd <- 2
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
          auth$n_free <- auth$n_free + n.auth.agents.remain
     } else {
          n.auth.agents.reserved <- auth$min_per_nbhd*n.nbhds.remain
          n.auth.agents.remain <- max( n.auth.agents.remain,  n.auth.agents.reserved)
          wts <- .HOSTNBHDSUM$n_agents[ .HOSTNBHDSUM$base_nbhd %in% remain_nbhds ]
          wts <- wts / sum(wts) 
          agent_partition <- sapply( wts, 
                                     FUN = function( x ){ 
                                          floor( x * (n.auth.agents.remain - n.auth.agents.reserved ) ) + auth$min_per_nbhd } 
                                     )
          names( agent_partition ) <- remain_nbhds
          # add whichever remain as free agents
          auth$n_free <- auth$n_free +  n.auth.agents.remain - sum( agent_partition )
          for( nbhd_name in remain_nbhds ){
               assigned_nbhds <- c( assigned_nbhds, rep(nbhd_name, times = agent_partition[nbhd_name]) )
               current_node <- c( current_node, 
                                  sample( AUTH_NBHD_NODES[[nbhd_name]], size = agent_partition[nbhd_name], replace = T) )
               }
     }
     
     assigned_nbhds <- c( assigned_nbhds, rep("free_", times = auth$n_free) )
     rem_nodes <- rownames( PLAYGROUND[!(PLAYGROUND$nbhd_name %in% assigned_nbhds), ] )
     if( length(rem_nodes) < auth$n_free ){
          rem_nodes <- rownames(PLAYGROUND)
     }
     current_node <- c( current_node, sample( rem_nodes, size = auth$n_free, replace = T) )
     
     df <- data.frame( assigned_nbhd = assigned_nbhds, current_node = current_node, 
                       stringsAsFactors = F )
     
     assign('AUTHORITIES_DF', df, .GlobalEnv)
     assign('AUTHORITIES', auth, .GlobalEnv)
     
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
     if( exists('AUTHORITIES_DF') ){
          if( nrow(AUTHORITIES_DF) > 0 ){
               auth_nbrs <- sapply( AUTHORITIES_DF$current_node, FUN = NbrIndices )
               auth_nbrs <- unique( do.call(c,auth_nbrs) )
          }
     }
     # collect hostile agent info who are near an authority agent
     auth_witnesses <- HOSTILES_DF[ HOSTILES_DF$current_node %in% auth_nbrs, ]
     witness_id <- auth_witnesses$uid
     witness_node <- auth_witnesses$current_node
     witness_nbhd <- auth_witnesses$base_nbhd
     witness_fam <- auth_witnesses$family
     witness_desig <- auth_witnesses$designation
     ## store these in the authority encounter log
     auth_sighting <- data.frame(  witness_id = witness_id, 
                                   witness_node = witness_node,
                                   witness_nbhd = witness_nbhd,
                                   witness_fam = witness_fam,
                                   witness_desig = witness_desig,
                                   stringsAsFactors = F )
     
     ## initialize the vectors that we want to remember
     ### for HOST_SIGHTING_HIST
     witness_id = character()
     witness_node = character()
     witness_nbhd = character()
     witness_fam = character()
     witness_desig = character()
     observed_id = character()
     observed_node = character()
     observed_nbhd = character()
     observed_fam = character()
     observed_desig = character()
     ### for HOSTILITY_HIST
     perp_id = character()
     perp_node = character()
     perp_nbhd = character()
     perp_fam = character()
     perp_desig = character()
     vict_id = character()
     vict_node = character()
     vict_nbhd = character()
     vict_fam = character()
     vict_desig = character()
     
     who_moved <- HOSTILES_DF$current_nbhd_desig != HOSTILES_DF$designation
     moved_agents <- HOSTILES_DF[ who_moved, ]
     unmoved_agents <- HOSTILES_DF[ !who_moved, ]
     
     n.moved <- nrow( moved_agents )
     for( mv_n in seq_maker(n.moved) ){
          agent <- moved_agents[mv_n, ]
          agent_nbrs <- as.character( agent$current_node )
          who_enemies <- moved_agents$current_node %in% agent_nbrs
          who_enemies <- who_enemies & (moved_agents$designation != agent$designation)
          moved_enemies <- moved_agents[ who_enemies, ] 
          n.moved_enemies <- nrow( moved_enemies )
          for( mv_n2 in seq_maker(n.moved_enemies) ){
               agent2 <- moved_enemies[mv_n2, ]
               # update the hostility encounter log
               witness_id = c( witness_id, agent$uid)
               witness_node = c( witness_node, agent$current_node )
               witness_nbhd = c( witness_nbhd, agent$base_nbhd )
               witness_fam = c( witness_fam, agent$family )
               witness_desig = c( witness_desig, agent$designation )
               observed_id = c( observed_id, agent2$uid )
               observed_node = c( observed_node, agent2$current_node )
               observed_nbhd = c( observed_nbhd, agent2$base_nbhd )
               observed_fam = c( observed_fam, agent2$family )
               observed_desig = c( observed_desig, agent2$designation )
               if( SHOOT(perp_agent = agent, vict_agent = agent2) ){
                    # update the hostility log
                    perp_id = c( perp_id, agent$uid )
                    perp_node = c( perp_node, agent$current_node )
                    perp_nbhd = c( perp_nbhd, agent$base_nbhd )
                    perp_fam = c( perp_fam, agent$family )
                    perp_desig = c( perp_desig, agent$designation )
                    vict_id = c( vict_id, agent2$uid )
                    vict_node = c( vict_node, agent2$current_node )
                    vict_nbhd = c( vict_nbhd, agent2$base_nbhd )
                    vict_fam = c( vict_fam, agent2$family )
                    vict_desig = c( vict_desig, agent2$designation )
               }
          }
          
          who_enemies <- unmoved_agents$current_node %in% agent_nbrs
          who_enemies <- who_enemies & (unmoved_agents$designation != agent$designation)
          unmoved_enemies <- moved_agents[ who_enemies, ] 
          n.unmoved_enemies <- nrow( unmoved_enemies )
          for( mv_n2 in seq_maker(n.unmoved_enemies) ){
               agent2 <- unmoved_enemies[mv_n2, ]
               # update the hostility encounter log agent -> agent2
               witness_id = c( witness_id, agent$uid)
               witness_node = c( witness_node, agent$current_node )
               witness_nbhd = c( witness_nbhd, agent$base_nbhd )
               witness_fam = c( witness_fam, agent$family )
               witness_desig = c( witness_desig, agent$designation )
               observed_id = c( observed_id, agent2$uid )
               observed_node = c( observed_node, agent2$current_node )
               observed_nbhd = c( observed_nbhd, agent2$base_nbhd )
               observed_fam = c( observed_fam, agent2$family )
               observed_desig = c( observed_desig, agent2$designation )
               # update the hostility encounter log agent2 -> agent
               witness_id = c( witness_id, agent2$uid)
               witness_node = c( witness_node, agent2$current_node )
               witness_nbhd = c( witness_nbhd, agent2$base_nbhd )
               witness_fam = c( witness_fam, agent2$family )
               witness_desig = c( witness_desig, agent2$designation )
               observed_id = c( observed_id, agent$uid )
               observed_node = c( observed_node, agent$current_node )
               observed_nbhd = c( observed_nbhd, agent$base_nbhd )
               observed_fam = c( observed_fam, agent$family )
               observed_desig = c( observed_desig, agent$designation )
               if( SHOOT(perp_agent = agent, vict_agent = agent2) ){
                    # update the hostility log agent -> agent2
                    perp_id = c( perp_id, agent$uid )
                    perp_node = c( perp_node, agent$current_node )
                    perp_nbhd = c( perp_nbhd, agent$base_nbhd )
                    perp_fam = c( perp_fam, agent$family )
                    perp_desig = c( perp_desig, agent$designation )
                    vict_id = c( vict_id, agent2$uid )
                    vict_node = c( vict_node, agent2$current_node )
                    vict_nbhd = c( vict_nbhd, agent2$base_nbhd )
                    vict_fam = c( vict_fam, agent2$family )
                    vict_desig = c( vict_desig, agent2$designation )
               }
               if( SHOOT(perp_agent = agent2, vict_agent = agent) ){
                    # update the hostility log agent2 -> agent
                    perp_id = c( perp_id, agent2$uid )
                    perp_node = c( perp_node, agent2$current_node )
                    perp_nbhd = c( perp_nbhd, agent2$base_nbhd )
                    perp_fam = c( perp_fam, agent2$family )
                    perp_desig = c( perp_desig, agent2$designation )
                    vict_id = c( vict_id, agent$uid )
                    vict_node = c( vict_node, agent$current_node )
                    vict_nbhd = c( vict_nbhd, agent$base_nbhd )
                    vict_fam = c( vict_fam, agent$family )
                    vict_desig = c( vict_desig, agent$designation )
               }
          }
     }
     
     host_sighting <- data.frame( witness_id = witness_id, 
                                  witness_node = witness_node,
                                  witness_nbhd = witness_nbhd,
                                  witness_fam = witness_fam,
                                  witness_desig = witness_desig,
                                  observed_id = observed_id,
                                  observed_node = observed_node,
                                  observed_nbhd = observed_nbhd,
                                  observed_fam = observed_fam,
                                  observed_desig = observed_desig )
     
     
     
     host_hist <- data.frame( perp_id = perp_id,
                              perp_node = perp_node,
                              perp_nbhd = perp_nbhd,
                              perp_fam = perp_fam, 
                              perp_desig = perp_desig,
                              vict_id = vict_id,
                              vict_node = vict_node,
                              vict_nbhd = vict_nbhd,
                              vict_fam = vict_fam,
                              vict_desig = vict_desig )
     
     
     
     if( nrow( host_hist ) > 0 ){
          # combine events where perps from the same nbhd, standing at the same node shoot at victs from the same nbhd
          host_hist <- dplyr::filter( .data = dplyr::group_by( .data = host_hist, perp_node, perp_nbhd, vict_nbhd ),
                                      row_number() == sample(seq_maker(n()), size = 1) )
          # combind events where multiple perps from the same nbhd shoot at the victims in the same node from the same nbhd
          host_hist <- dplyr::filter( .data = dplyr::group_by( .data = host_hist, vict_node, vict_nbhd, perp_nbhd ),
                                      row_number() == sample(seq_maker(n()), size = 1) )
          host_hist <- as.data.frame( lapply( host_hist, FUN = as.character ), stringsAsFactors = F )
     }
     ifelse( nrow(host_hist) > 0,
             host_hist$timestep <- as.integer(timestep),
             host_hist$timestep <- integer() )

     
     if( nrow( host_sighting ) > 0) {
          # combine sightings of agents from the same nbhd at the same node witnessing enemies from the same family
          host_sighting <- dplyr::filter( .data = dplyr::group_by( .data = host_sighting, witness_node, witness_nbhd, observed_fam ),
                                          row_number() == sample(seq_maker(n()), size = 1) )
          host_sighting <- as.data.frame( lapply( host_sighting, FUN = as.character ), stringsAsFactors = F  )
     }
     ifelse( nrow(host_sighting) > 0,
             host_sighting$timestep <- as.integer(timestep),
             host_sighting$timestep <- integer() )
     
     if( nrow( auth_sighting ) > 0 ){
          # combine sightings of agents from the same nbhd at the same node witnessing authorities
          auth_sighting <- dplyr::filter( .data = dplyr::group_by( .data = auth_sighting, witness_node, witness_nbhd ),
                                          row_number() == sample(seq_maker(n()), size = 1) )
          auth_sighting <- as.data.frame( lapply( auth_sighting, FUN = as.character ), stringsAsFactors = F )
     }
     ifelse( nrow(auth_sighting) > 0,
             auth_sighting$timestep <- as.integer(timestep),
             auth_sighting$timestep <- integer() )
     
     if( nrow(auth_sighting) > 0 ){
          # check if the auth_sighting prevented a hostility
          if( AGGL_AUTHORITY ){
               # record which interaction prevented hostility
               auth_to_host <- auth_sighting$witness_node %in% host_hist$perp_node
               auth_to_host <- auth_to_host | (auth_sighting$witness_node %in% host_hist$vict_node)
               auth_sighting$prevented_host <- auth_to_host
               # remove those prevented hostilities from the record (they never happened)
               host_to_auth <- host_hist$perp_node %in% auth_sighting$witness_node
               host_to_auth <- host_to_auth | (host_hist$vict_node %in% auth_sighting$witness_node)
               host_hist <- host_hist[ - which( host_to_auth ), ]
          } else {
               # record which interaction prevented hostility
               auth_to_host <- auth_sighting$witness_node %in% host_hist$perp_node
               auth_sighting$prevented_host <- auth_to_host
               # remove those prevented hostilities from the record (they never happened)
               host_to_auth <- host_hist$perp_node %in% auth_sighting$witness_node
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
               host_sighting <- rbind( HOST_SIGHTING_DF, host_sighting, make.row.names = F )
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
     saveRDS( list(LATTICE = LATTICE, 
                   HOSTILES = HOSTILES, 
                   PLAYGROUND = PLAYGROUND, 
                   .PGSUM = .PGSUM, 
                   .ETHSUM = .ETHSUM), 
              file = file_path )
}

OpenPlayground <- function( file_path ){
     pg <- readRDS( file_path )
     LATTICE <<- pg$LATTICE
     HOSTILES <<- pg$HOSTILES
     PLAYGROUND <<- pg$PLAYGROUND
     .PGSUM <<- pg$.PGSUM
     .ETHSUM <<- pg$.ETHSUM
     cat(paste0("Succesfully loaded file.\n",
               "LATTICE, HOSTILES, PLAYGROUND, .PGSUM, .ETHSUM are now in your global environment.\n") )
}


#################################### VISUALIZATION ####################################
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

.VisHostiles <- function( alph = .5 ){
     if( !exists('HOSTILES_DF') ){
          cat('Must create HOSTILES_DF using CreateHostilesDF before visualization.\n')
          return( F )
     }
     hostplot <- PLAYGROUND[HOSTILES_DF$current_node, c('row','col','family')]
     max_row <- max(PLAYGROUND$row)
     hostplot <- dplyr::mutate(hostplot, plot_row = max_row + 1 - row )
     
     
     fams <- unique(hostplot$family)
     colourCount = length(fams)
     getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,"Reds"))
     hostplot$color <- character(nrow(hostplot))
     for( i in 1:colourCount ){
          hostplot$color[hostplot$family == fams[i]] <- getPalette(colourCount)[i]
     }
     
     g <-  ggplot2::geom_point( data = hostplot, 
                                   ggplot2::aes(x = col, y = plot_row),
                                   color = hostplot$color,
                                   alpha = alph )
     
     return( g ) 
}

.VisAuthorities <- function( alph = .7 ){
     if( !exists('AUTHORITIES_DF') ){
          cat('Must create AUTHORITIES_DF using CreateAuthoritiesDF before visualization.\n')
          return( F )
     }
     authplot <- PLAYGROUND[AUTHORITIES_DF$current_node, c('row','col','family')]
     max_row <- max(PLAYGROUND$row)
     authplot <- dplyr::mutate(authplot, plot_row = max_row + 1 - row )
     
     g <- ggplot2::geom_point( data = authplot, 
                                   ggplot2::aes(x = col, y = plot_row),
                                   color = "red",
                                   alpha = alph )
     
     return( g )
}

ViewPlayground <- function( alph = 1 ){
     g <- .VisPlayground( alph = alph ) + ggplot2::ggtitle( "PLAYGROUND" )
     return( g )
}

ViewHostiles <- function( alph = .25 ){
     g <- .VisPlayground( alph = 1 ) + .VisHostiles( alph = alph )
     return( g )
}

ViewAuthorities <- function( alph = .75 ){
     g <- .VisPlayground( alph = .33 ) + .VisAuthorities( alph = alph )
     return( g )
}

ViewAgents <- function(){
     g <- .VisPlayground( alph = .33 ) + .VisHostiles( alph = .5 ) + .VisAuthorities( alph = .75 )
     return( g )
}

#######################################################################################
#######################################################################################
############################### SECTION: TUNABLE FUNCTIONS ############################
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
          return( max(TENSION_MATRIX[grp1,grp2], TENSION_MATRIX[grp2,grp1]) )
     } else if( mode == "min" ){
          return( min(TENSION_MATRIX[grp1,grp2], TENSION_MATRIX[grp2,grp1]) )
     } else {
          return( TENSION_MATRIX[from_grp, towards_grp] )
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
     tens <- NBHD_TENSION( perp_nbhd, vict_nbhd )
     # find the probability of shooting
     p_shoot <- 1 - exp(- shoot_factor * (perp_viol * grp_tens)) #basically = lam * perp_viol * grp_tens when small
     # return decision 
     return( sample( c(TRUE,FALSE), size = 1, prob = c(p_shoot,1-p_shoot) ) )
}

# for ATTACK
AttackDecision <- function( hfam, allow_multiple_attacks = F ){
     # requires DOMINANCE MATRIX
     # returns a logical vector indexed by family names, T => attack
     will_attack_at <- 5
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

############################### SECT: TUNABLE FUNCTIONS ###############################
############################### SUBSECTION: PMFers ####################################
#######################################################################################

StdHostMovesNbhd <- function( hnbhd,
                              home_wt, 
                              ether_wt,
                              ether_nbhd_wts ){
     
     den <- home_wt + ether_wt
     if( ! ( home_wt < 0 | ether_wt < 0 | den == 0) ){
          homep <- home_wt / den
          etherp <- ether_wt / den
     } else {
          homep <- 1
          etherp <- 0
     }
     if( !( is.null(ether_nbhd_wts) | any(ether_nbhd_wts < 0) | sum(ether_nbhd_wts) == 0 ) ){
          ether_nbhd_ps <- ether_nbhd_wts / sum( ether_nbhd_wts )
     } else {
          ether_nbhd_ps <- c(ether = 1)
     }
     pmf <- c(homep, etherp * ether_nbhd_ps)
     names( pmf ) <- c( hnbhd, names(ether_nbhd_ps) )
     
     where <- sample( names(pmf), size = 1, prob = pmf )
     
     return( where )
     
}


StdHostMovesNode <- function( agent,
                              nbhd_to,
                              timestep,
                              ether_nbhd_names = 'ether' ){
     
     base_nbhd <- agent$base_nbhd
     
     ## if staying home, stay home
     if( nbhd_to == base_nbhd ) return( agent$base_node )
     
     ## how far into the nbhd
     dists <- unique( PLAYGROUND[PLAYGROUND$nbhd_name == nbhd_to, paste0('dist_from_', base_nbhd)] )
     dists_pmf <- DIST_PMFer( dists )
     # if there was an error, stay home
     if( is.null( dists_pmf ) ) return( agent$base_node )
     dist_in <- sample( dists, size = 1, prob = dists_pmf )
     
     ## which nodes are accessible at that distance
     red <- PLAYGROUND[ PLAYGROUND$nbhd_name == nbhd_to, ]
     red <- red[ red[[paste0('dist_from_',base_nbhd)]] == dist_in, ]
     nodes <- rownames( red )
     
     if( is.null( auth_pmf <- AVOID_AUTH_NODES_PMFer(timestep = timestep,
                                                      agent = agent,
                                                      nbhd_to = nbhd_to,
                                                      dist_in = dist_in,
                                                      nodes = nodes ) ) ){
          # if something went wrong, stay home
          n.nodes <- length(nodes)
          auth_pmf <- rep( 1, n.nodes ) / n.nodes
          names( auth_pmf ) <- nodes
     }
     
     if( nbhd_to %in% ether_nbhd_names ){
          if( is.null( enemy_pmf <- AVOID_HOST_NODES_PMFer(timestep = timestep,
                                                                   agent = agent,
                                                                   nbhd_to = nbhd_to,
                                                                   dist_in = dist_in,
                                                                   nodes = nodes ) ) ){
               # if something went wrong, stay home
               n.nodes <- length(nodes)
               enemy_pmf <- rep( 1, n.nodes ) / n.nodes
               names( enemy_pmf ) <- nodes
          }
          use_pmf <- auth_pmf * enemy_pmf
          use_pmf <- use_pmf / sum( use_pmf )
     } else {
          use_pmf <- auth_pmf
     }
     
     return( sample( nodes, size = 1, prob = use_pmf ) )
}

StdHostAgentMover <- function( agent,
                               timestep,
                               home_wt = 8,
                               ether_wt = 1,
                               ether_nbhd_wts = NULL){
     moved_agent <- agent
     hnbd <- moved_agent$base_nbhd
     moved_agent$current_nbhd <- StdHostMovesNbhd( hnbd, home_wt, ether_wt, ether_nbhd_wts )
     moved_agent$current_node <- StdHostMovesNode( moved_agent, moved_agent$current_nbhd, timestep )
     return( moved_agent )
}


# for DIST_PMFer
PMFer_ByDists <- function( dists,
                           dfun = function( x ){ x^2 },
                           offset = 0.5 ){
     dists_pmf <- 1 / (dfun(dists) + offset)
     if( (den <- sum(dists_pmf)) == 0 ) return( NULL )
     dists_pmf <- dists_pmf / den
     names( dists_pmf ) <- as.character( dists )
     return( dists_pmf )
}

# for AVOID_AUTH_NODES_PMFer
PMFer_HideFromAuthNodes <- function( agent, 
                                     nbhd_to,
                                     dist_in,
                                     nodes,
                                     timestep,
                                     wt_offset = 1,
                                     mem_exp_fac = 2){
     
     
     
     base_nbhd <- agent$base_nbhd
     
     CalculateAuthSightingNodeWts <- function( ){
          # requires AUTH_SIGHTING_HIST, PLAYGROUND
          if( ! exists( 'AUTH_SIGHTING_HIST' ) ) return(NULL)
          most_recent <- max( AUTH_SIGHTING_HIST$timestep )
          wts <- NULL
          auth <- AUTH_SIGHTING_HIST[AUTH_SIGHTING_HIST$witness_nbhd == base_nbhd, ]
          for( timestep in unique(auth$timestep) ){
               cent_nodes <- auth[auth$timestep == timestep, ]$witness_node
               auth_sighting_nodes <- NULL
               for( node_index  in cent_nodes ){
                    auth_sighting_nodes <- c( auth_sighting_nodes, NbrIndices( node_index ) )
               }
               print(auth_sighting_nodes)
               print(nodes)
               auth_sighting_nodes <- intersect( auth_sighting_nodes, nodes )
               for( node_index in auth_sighting_nodes ){
                    if( node_index %in% names(wts) ){
                         wts[node_index] <- wts[node_index] + mem_exp_fac^{-(most_recent - timestep)}
                    } else {
                         wts <- c(wts, mem_exp_fac^{-(most_recent - timestep)})
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
PMFer_HideFromEnemyNodesByDist <- function( agent,
                                            nbhd_to,
                                            dist_in,
                                            nodes,
                                            timestep ){
     
     # currently only requires .ETHSUM
     run_anew <- function(){
          dists <- .ETHSUM[nodes, -which( names(.ETHSUM) %in% paste0('dist_from_', agent$designation) ) ]
          names( dists ) <- nodes
          dists_pmf <- dists + 10
          if( (den <- sum(dists_pmf)) == 0 ) return( NULL )
          return( dists_pmf / den )
     }
     
     ## use the timestep to store hashing
     base_nbhd <- agent$base_nbhd
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
                                                mem_exp_fac = 2 ){
     
     
     # used in PMFer_HideFromEnemyNodesBySightings isn't based off sightings
     CalculateHostSightingNodeWts <- function( ){
          # requires HOST_SIGHTING_HIST, PLAYGROUND
          if( ! exists( 'HOST_SIGHTING_HIST' ) ) return(NULL)
          most_recent <- max( HOST_SIGHTING_HIST$timestep )
          wts <- NULL
          auth <- HOST_SIGHTING_HIST[HOST_SIGHTING_HIST$witness_nbhd == base_nbhd, ]
          for( timestep in unique(auth$timestep) ){
               cent_nodes <- auth[auth$timestep == timestep, ]$witness_node
               enemy_sighting_nodes <- NULL
               for( node_index  in cent_nodes ){
                    enemy_sighting_nodes <- c( enemy_sighting_nodes, NbrIndices( node_index ) )
               }
               enemy_sighting_nodes <- intersect( enemy_sighting_nodes, nodes )
               for( node_index in enemy_sighting_nodes ){
                    if( node_index %in% names(wts) ){
                         wts[node_index] <- wts[node_index] + mem_exp_fac^{-(most_recent - timestep)}
                    } else {
                         wts <- c(wts, mem_exp_fac^{-(most_recent - timestep)})
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


# for NBHD_ATTACK_PMFer
PMFer_NbhdToAttack <- function( attacker_fam, 
                                attackee_fam,
                                mem_exp_fac ){
     
     # requires HOSTILITY_HIST and PLAYGROUND
     
     # nbhds belonging to the family who will be the victim of the attack
     attackee_nbhds <- unique( HOSTILES_DF$base_nbhd[ HOSTILES_DF$family == attackee_fam ] )
     print( attackee_nbhds )
     CalculateAttackNbhdWts <- function( ){
          # requires HOSTILITY_HIST, PLAYGROUND
          if( ! exists( 'HOSTILITY_HIST' ) ) return(NULL)
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
                         wts[nb] <- wts[nb] + mem_exp_fac^{-(most_recent - timestep)}
                    } else {
                         wts <- c(wts, mem_exp_fac^{-(most_recent - timestep)})
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
PMFer_AssignedAgent <- function( agent, timestep, wt_offset, mem_exp_fac ){
     # requires AUTH_NBHD_NODES
     CalculateHotNodeWts <- function( ){
          # requires HOSTILITY_HIST, PLAYGROUND
          if( ! exists( 'HOSTILITY_HIST' ) ) return(NULL)
          most_recent <- max( HOSTILITY_HIST$timestep )
          wts <- NULL
          for( timestep in unique(HOSTILITY_HIST$timestep) ){
               cent_nodes <- HOSTILITY_HIST[HOSTILITY_HIST$timestep == timestep, ]$perp_node
               hot_nodes <- NULL
               for( node_index  in cent_nodes ){
                    hot_nodes <- c( hot_nodes, NbrIndices( node_index ) )
               }
               for( node_index in hot_nodes ){
                    if( node_index %in% names(wts) ){
                         wts[node_index] <- wts[node_index] + mem_exp_fac^{-(most_recent - timestep)}
                    } else {
                         wts <- c(wts, mem_exp_fac^{-(most_recent - timestep)})
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


#################################### TUNABLE VARIABLES #################################
LATTICE <- RectLatticeMaker( nrows = 350, ncols = 150 ) 

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
     n_assigned = 100,
     n_free = 10,
     assignment_radius = 5,
     min_per_nbhd = 2
)

### Any other type of nbhd
OTHER_NBHDS <- list()

### Number of authority agents
N_AUTHORITIES <- 100

### Nearby neighbor radius
NBR_RADIUS <- 1

### Use "walled neighborhoods" during agent interaction
WALLED_NBHDS <- T

### Agglomerative authority interaction
AGGL_AUTHORITY <- T

### Run in VERBOSE mode
VERBOSE <- T

##########################################
## GLOBAL VARIABLES FROM TUNABLE FUNCTIONS
##########################################

# METRIC must accept node1 and node2 as entries
METRIC <- function(node1, node2, ...){ TaxiNodeNodeDist(node1, node2 ) }

# NBHD_TENSION must accept from_nbhd and towards_nbhd as entries
NBHD_TENSION <- function(from_nbhd, towards_nbdh, ...){ TensionCalculator(from_nbhd, towards_nbhd, mode = 'max') }

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

# ATTACK must accept hfam; returns chr of family to attack, or NULL if none
ATTACK <- function( hfam, ... ){ AttackDecision( hfam, allow_multiple_attacks = F ) }

### PMFers
# ASSIGNED_AGENT_PMFer must accept agent and timestamp
ASSIGNED_AGENT_PMFer <- function( agent, timestep, ... ){ PMFer_AssignedAgent( agent, 
                                                                               timestep, 
                                                                               wt_offset = 1, 
                                                                               mem_exp_fac = 2 ) }
# NBHD_ATTACK_PMFer must accept attacker_fam and attackee_fam
NBHD_ATTACK_PMFer <- function( attacker_fam, attackee_fam, ... ){ PMFer_NbhdToAttack( attacker_fam, 
                                                                                      attackee_fam,
                                                                                      mem_exp_fac = 1.5 ) }

# DIST_PMFer must accept dists
DIST_PMFer <- function( dists, ... ){ PMFer_ByDists(dists, 
                                                    dfun = function(x){ x^2},
                                                    offset = 0.5 ) }

# HIDE_FROM_ENEMY_NODES_PMFer
AVOID_ENEMY_NODES_PMFer <- function( agent, nbhd_to, dist_in, nodes, timestep, ...){ 
     PMFer_HideFromEnemyNodesByDist( agent,
                                     nbhd_to,
                                     dist_in,
                                     nodes,
                                     timestep )
}

# AVOID_AUTH_NODES_PMFer
AVOID_AUTH_NODES_PMFer <- function( agent, nbhd_to, dist_in, nodes, timestep, ...){ 
     PMFer_HideFromAuthNodes( agent, 
                              nbhd_to,
                              dist_in,
                              nodes,
                              timestep,
                              wt_offset = 1,
                              mem_exp_fac = 2)
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
     if( !exists('LATTICE') ){
          cat(paste("ERROR: LATTICE has not been defined.","There is no default list for LATTICE\n", sep = " "))
          return( F )
     }
     if( !exists('HOSTILES') ){
          cat(paste("ERROR: HOSTILES list not found.","There is no default list for HOSTILES.\n", sep = " "))
          return( F )
     }
     if( !exists('AUTHORITIES') ){
          cat("- AUTHORITIES list not found. Setting to default list with 275 agents.\n")
          l <- list( n_assigned = 250, n_free = 25, assign_nbhd_by_population = T )
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
          cat("Error creating HOSTILES_DF dataframe via CreateHostilesDF.\nSetup Failed.")
          return(F)
     }
     if( !CreateHostNbhdSummary() ){ 
          cat("Error creating .HOSTNBDHSUM dataframe via CreateHostNbhdSummary.\nSetup Failed.")
          return(F)
     }
     if(!CreateAuthNbhdNodes()){
          cat(paste0("Error creating AUTH_NBHD_NODES, necessasary for authority agents in model.\n",
                     "No authority agents will be added to model.\n") )
          CreateEmptyAuthoritiesDF()
     } else if(!CreateAuthoritiesDF()){ 
          cat(paste0("Error creating AUTHORITIES_DF, necessasary for authority agents in model.\n",
                     "No authority agents will be added to model.\n") )
          CreateEmptyAuthoritiesDF()
     }
     return(T)
}
