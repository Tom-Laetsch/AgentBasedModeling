require(dplyr)
require(ggplot2)

"
We construct an Agent Based Model (ABM), partitioning the agents into N+1
groups: N hostile groups; 1 police group. At each timestep (tick), each 
agent will move according to a probability distribution on the ABM lattice
constructed below. When agents from different hostile groups are within a 
certain radius of each other, they will decide whether or not to act hostile
towards each other, which roughly can be described as: if as police agent is 
near, no hostility will take place; however, if no police agent is near then 
the opposing agents will sample from a Bernoulli distribution (with distribution
dependent on previous hostilities between these groups) and decide whether to
enact a hostility on the opposing group based on the sampling. 
"

################################################################################
############################ GRID LATTICE CONSTRUCTION #########################
"
In this section we house a few useful functions for creating the grid lattice on 
which the agents list, and defining metrics on this lattice, which is where are
used in determination of neighbors to a given agent and potentially affecting the
movement distribution of each agent. 

Some conventions and notes:
- Each node on the lattice will be referred to as a 'block'
- The lattice is conceptualized as a two-column matrix where adjacent blocks (nodes)
     are then realized as adjacent entries within this matrix. For example, 
     with a very simple square grid having only four blocks, the grid would be 
     conceptualized as a square 2x2 matrix.
- Each hostile group is given a neighborhood on the lattice, which is stored as 
     a collection of row/col numbers referring to the row/column that neighborhood
     inhabits of the grid matrix. 
- Based on a tunable radius (see NEARBY_RADIUS), two agents on the grid will 
     be considered interacting neighbors if the distance (see BlcBlcDist) between
     them is below this radius. Brute force calculations of nearest neighbors is 
     computationally expensive, so to speed this up, when the grid is created and
     a nearby radius determined, we create the list NEARBY_LIST, with each block 
     corresponding to a list entry which containing the adjacent blocks within the 
     radius. While this uses more memory, it dramatically speeds up the calcuations 
     (e.g., the brute force HostilityUpdater_old is ~60x slower than this newer 
     technique used in HostilityUpdater). 
"

RectNbhdMaker <- function( tl_br ){
"
RectNbhdMaker: Creates a matrix of coordinates of a rectangular region specified by inputs. 

- input: tl_br = a numeric 4-vector where the first entry is the top-left row number,
     the second is the top-left col number, the third is the bottom-right row number, 
     and the fourth is the bottom-right column number. tl_br = c( tlr, tlc, brr, brc )

- output: A two-column matrix with the first column being the row-numbers and the second
     being the corresponding column numbers for the lattice region contained within the 
     passed top-left & bottom-right numbers.
"
        tlr = tl_br[1]
        tlc = tl_br[2]
        brr = tl_br[3]
        brc = tl_br[4]
        rlen <- brr - tlr + 1
        clen <- brc - tlc + 1
        nbh <- matrix( c(rep( tlr:brr, each=clen ), rep(tlc:brc, times = rlen)), ncol = 2 )
        colnames(nbh) <- c("row","col")
        return( nbh )
}

GridMaker <- function( nrows, ncols ){
"
GridMaker: Creates a rectangular lattice matrix starting at row=1, col=1 with dimensions specified by inputs.

- input: two integer values nrows, ncols which specify the number of rows and columns, resp., to use
     when making the lattice.

- output: A two-column matrix with the first column being the row-numbers and the second being the 
     corresponding column-numbers for the lattice grid starting at top-left: (row,col)=(1,1) and 
     convering the lattice points with bottom-right: (row,col)=(nrow,ncol). 
"
     return( RectNbhdMaker( c(topl_row = 1, topl_col = 1, botr_row = nrows, botr_col = ncols)) )
}



######### GRID METRICS 
# distance between blocks
BlcBlcDist <- function( b1, b2 ){
        return( max( abs(b1[1]-b2[1]), abs(b1[2] - b2[2]) ) )
}
# distance between a block and a nbhd
BlcNbhdDist <- function( b, nbhd_coords ){
        return( min( apply( X=nbhd_coords, FUN=function(x){ BlcBlcDist( x,b ) }, MARGIN=1 ) ) )
}
# distance between two nbhds
NbhdNbhdDist <- function( nbhd_coords1, nbhd_coords2 ){
        return( min( apply( nbhd_coords1, MARGIN=1, FUN=function(x){ BlcNbhdDist( x, nbhd_coords2 ) } ) ) )
}
#tests if two blocks are neighbors 
AreNeighbors <- function( b1, b2, radius ){
     return( BlcBlcDist(b1,b2) <= radius )
}

# create a list containing areas of the grid which we search for neighbors 
NearbyListNameMaker <- function( blc ){ 
     #helper function to give consistent list names in ConstructNearbyList 
     #and elsewhere when we need to grab elements from the list
     return( paste0("r",blc[1],"_c",blc[2]) ) 
}
ConstructNearbyList <- function( grid_layout, nearby_radius ){
     nearby_list <- list()
     for( i in 1:nrow(grid_layout) ){
          bi <- as.matrix( grid_layout[ i, c('row','col') ] )
          ci <- NearbyListNameMaker(bi)
          nearby_list[[ci]] <- list( centroid = bi, nearby = bi )
     }
     nombres <- names(nearby_list)
     n_noms <- length(nombres)
     for( i in 1:(n_noms-1) ){
          ci <- nombres[ i ]
          bi <- nearby_list[[ci]]$centroid
          for( j in (i+1):n_noms ){
               cj <- nombres[ j ]
               bj <- nearby_list[[cj]]$centroid
               if( AreNeighbors(bi,bj,nearby_radius) ){
                    nearby_list[[ci]]$nearby <- rbind(nearby_list[[ci]]$nearby, bj)
                    cj <- NearbyListNameMaker( bj )
                    nearby_list[[cj]]$nearby <- rbind(nearby_list[[cj]]$nearby, bi)
               }
          }
     }
     return(nearby_list)
}


########### INITIALIZATION STUFF
SymMatrixMaker <- function( names, init_val ){
        n <- length(names)
        mat <- matrix( rep( init_val, times = n*n ), nrow = n, ncol = n )
        colnames( mat ) <- names
        rownames( mat ) <- names
        return(mat)
}

AngerMatrixInit <- function( names, min_anger ){
        mat <- SymMatrixMaker( names, min_anger )
        #zeros down the diagonal (no anger within same group)
        mat[ row(mat) == col(mat) ] <- 0
        return(mat)
}

HostilityMatrixInit <- function( names ){
        # initialize so that every group has no hostilities against each other
        return( SymMatrixMaker(names, init_val = 0) )
}

AreaPMFInit <- function( group_names,
                         area_names,
                         grid_layout,
                         within_grp_p = 0.7, 
                         opposing_grp_p = 0.1,
                         ether_p = 0.2){
        area_pmf <- SymMatrixMaker( area_names, init_val = ether_p )
        for( from_grp in group_names ){
                movement_factor <- numeric()
                for( to_grp in group_names ){
                        if( from_grp == to_grp ){
                                area_pmf[from_grp,to_grp] <- within_grp_p
                        }
                        else{
                             from_nbhd <- as.matrix(grid_layout[grid_layout$area == from_grp, ][,c('row','col')])
                             to_nbhd <- as.matrix(grid_layout[grid_layout$area == to_grp, ][,c('row','col')])
                             movement_factor[to_grp] <- 1.0 / NbhdNbhdDist( from_nbhd, to_nbhd )
                        }
                }
                movement_factor <- opposing_grp_p * movement_factor / sum(movement_factor)
                area_pmf[from_grp, names(movement_factor)] <- movement_factor
        }
        return(area_pmf)
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

###### GROUP NEIGHBOR INTERACTION STUFF
AngerCalculator <- function( n_hostilities, prev_anger, min_anger ){
     
     if( n_hostilities == 0 ){
          forget_fac <- 0.9
          prev_fac <- prev_anger * forget_fac
          return( max( prev_fac, min_anger ) )
     }
     else if( n_hostilities == 1 ){
          hostility_fac <- 0.5
     }
     else if( n_hostilities == 2){
          hostility_fac <- 0.75
     }
     else if( n_hostilities >= 3){
          hostility_fac <- 1.0
     }
     
     return( hostility_fac + (1 - hostility_fac) * prev_anger )
     
}

AngerMatrixUpdater <- function( anger_matrix, hostility_matrix, min_anger ){
        n <- nrow(hostility_matrix)
        for( i in 1:n ){
                for( j in 1:n ){
                     if( i!=j ){
                          # note j,i vs i,j below: that way host&tens matrices are aligned the same: row=to, col=from
                          anger_matrix[i,j] <- AngerCalculator( hostility_matrix[j,i], anger_matrix[i,j], min_anger )
                     }
                }
        }
        return( anger_matrix )
}


TensionCalculator <- function( grp1, grp2, anger_matrix ){
        # given the asymmetric anger_matrix, return the symmetric tension between groups
     return( anger_matrix[grp1,grp2])
        return( max(anger_matrix[grp1,grp2], anger_matrix[grp2,grp1]) )
}

ShootDecision <- function( perp_grp, perp_viol, vict_grp, anger_matrix, lam ){
     grp_tens <- TensionCalculator( perp_grp, vict_grp, anger_matrix )
     p_shoot <- 1 - exp(- lam * (perp_viol * grp_tens))
     return( sample( c(TRUE,FALSE), size = 1, prob = c(p_shoot,1-p_shoot) ) )
}

HostilesNotNearAuthority <- function( authority_df, hostiles_df, nearby_radius ){
     ret_rows <- rep(TRUE, times = nrow(hostiles_df) )
     for( officer_n in seq_len(nrow(authority_df)) ){
          officer_b <- authority_df[officer_n, c('row','col')]
          for( hostile_n in seq_len(nrow(hostiles_df)) ){
               hostile_b <- hostiles_df[hostile_n, c('row','col') ]
               if( BlcBlcDist(officer_b, hostile_b) <= nearby_radius ){
                    ret_rows[hostile_n] <- FALSE
               }
          }
     }
     return( hostiles_df[ret_rows, ] )
}

#_old implies outdated and not used
HostilityUpdater_old <- function( hostility_record_df,
                              hostiles_df, 
                              authority_df,
                              anger_matrix,
                              current_timestep,
                              nearby_radius,
                              lam ){
     
     # step 1: remove consideration of hostile agents near police
     use_hostiles <- HostilesNotNearAuthority( authority_df, hostiles_df, nearby_radius )
     # step 2:
     rem_groups <- unique( as.character(use_hostiles$group) )
     # step 3:
     if( length(rem_groups) <= 1 ){ return(hostility_record_df) }
     for( grp in rem_groups ){
          perps_l <- use_hostiles$group == grp
          perps <- use_hostiles[ perps_l, ]
          victs <- use_hostiles[ !perps_l, ]
          for( perp_n in seq_len(nrow(perps)) ){
               perp <- perps[ perp_n, ]
               for( vict_n in seq_len(nrow(victs)) ){
                    vict <- victs[vict_n, ]
                    # check if perp and vict are neighbors
                    if( BlcBlcDist(b1 = perp[c('row','col')], 
                                   b2 = vict[c('row','col')]) <= nearby_radius ){
                         # if so, decide if perp shoots at vict
                         perp_grp <- as.character(perp$group)
                         vict_grp <- as.character(vict$group)
                         if( ShootDecision(perp_grp = perp_grp, 
                                           perp_viol = perp$violence_fac, 
                                           vict_grp = vict_grp, 
                                           anger_matrix = anger_matrix,
                                           lam = lam) ){
                              new_row <- data.frame( time = current_timestep, 
                                                     perp_grp = perp_grp,
                                                     vict_grp = vict_grp,
                                                     row = perp$row,
                                                     col = perp$col )
                              hostility_record_df <- rbind( hostility_record_df, new_row )
                         }
                    }
               }
          }
     }
     
     return( hostility_record_df )
}


HostilityUpdater <- function( current_timestep,
                              hostility_record_df, 
                              hostiles_df, 
                              authority_df, 
                              nearby_list,
                              anger_matrix,
                              lam ){
     agents_nearby <- function( agents_rc_mat, nearby_mat ){
          return( apply(agents_rc_mat, MARGIN = 1, 
                        FUN = function(row1){ 
                                             any( apply(nearby_mat, 1, function(row2){ all(row1 == row2) }) ) 
                                            }
                       ) 
                )
     }
     
     hostiles_mat <- as.matrix( hostiles_df[c('row','col')] )
     
     to_ignore_cent_names <- apply( as.matrix( unique( authority_df[ c('row','col')] ) ), 
                                    MARGIN = 1,
                                    FUN = NearbyListNameMaker )
     to_ignore <- NULL
     for( cent in to_ignore_cent_names ){
          to_ignore <- rbind( to_ignore, nearby_list[[cent]]$nearby )
     }
     
     if(!is.null(to_ignore)){
          to_ignore_names <- apply( unique( to_ignore ), 
                                    MARGIN = 1,
                                    FUN = NearbyListNameMaker )
          
     }else{ to_ignore_names <- NULL }
     to_use_cent_names <- apply( unique( hostiles_mat ) , 
                                 MARGIN = 1,
                                 FUN = NearbyListNameMaker )
     to_use_cent_names <- setdiff( to_use_cent_names, to_ignore_names )
     for( cent in to_use_cent_names ){
          nls <- nearby_list[[cent]]
          cent_who <- apply( hostiles_mat, MARGIN = 1, FUN = function( row ){ all( row == nls$centroid ) } )
          near_who <- agents_nearby( hostiles_mat, nls$nearby )
          cent_hosts <- hostiles_df[ cent_who, ]
          near_hosts <- hostiles_df[ near_who, ]
          cent_grps <- unique( cent_hosts$group )
          near_grps <- unique( near_hosts$group )
          if( length(near_grps) <= 1 ) next
          for( perp_grp in cent_grps ){
               perps <- cent_hosts[ cent_hosts$group == perp_grp, ]
               victs <- near_hosts[ near_hosts$group != perp_grp, ]
               for( perp_n in seq_len(nrow(perps)) ){
                    perp <- perps[perp_n, ]
                    for( vict_n in seq_len(nrow(victs)) ){
                         vict <- victs[vict_n, ]
                         if( ShootDecision(perp_grp = perp$group, 
                                           perp_viol = perp$violence_fac, 
                                           vict_grp = vict$group, 
                                           anger_matrix = anger_matrix,
                                           lam = lam) ){
                              new_row <- data.frame( time = current_timestep, 
                                                     perp_grp = perp$group,
                                                     vict_grp = vict$group,
                                                     row = perp$row,
                                                     col = perp$col )
                              hostility_record_df <- rbind( hostility_record_df, new_row )
                         }
                    }
               }
          }
     }
     return( hostility_record_df )
}

HostilityMatrixMaker <- function( group_names, hostility_record_df, timestep ){
     hostility_matrix <- SymMatrixMaker(names = group_names, init_val = 0)
     hostilities <- hostility_record_df[ hostility_record_df$time == timestep, ] %>%
                    group_by( perp_grp, vict_grp ) %>%
                    summarise( n = n() )
     for( rown in seq_len(nrow(hostilities)) ){
          pgrp <- as.character(hostilities[rown, ]$perp_grp)
          vgrp <- as.character(hostilities[rown, ]$vict_grp)
          n <- as.numeric(hostilities[rown, ]$n)
          hostility_matrix[ pgrp, vgrp ] <- n
     }
     return(hostility_matrix)
}


############### MOVEMENT FUNCTIONS

SingleAgentMoveUpdater <- function( agent, area_names, grid_layout ){
     # step 1: choose which area to move to based on pmf
     area <- sample( area_names, size = 1, prob = agent[ area_names ])
     # step 2: choose where inside that area to move
     # currently set to uniform, but can change
     idxs <- as.numeric(rownames(grid_layout[grid_layout$area == area, ]))
     idx <- sample( idxs, size = 1 )
     # step 3: return the location
     return( c( row = grid_layout[idx,'row'], col = grid_layout[idx, 'col'] ) )
}

MoveUpdater <- function( agent_df, area_names, grid_layout ){
     for( i in 1:nrow(agent_df) ){
          move <- SingleAgentMoveUpdater( agent_df[i, ], area_names, grid_layout )
          agent_df[i, 'row'] <- move['row']
          agent_df[i, 'col'] <- move['col']
     }
     return( agent_df )
}

############### CONSTRUCTION 

GRID <- GridMaker( nrows = 12, ncols = 12 )
     #RectNbhdMaker(c(topl_row = 1, topl_col = 1, botr_row = 12, botr_col = 12))

MIN_ANGER = 0.1 # minimal anger between hostile groups
NEARBY_RADIUS <- 1 # = n => check current block + all blocks within distance of n

HOSTILES <- list(
     hgrp1 = list( nbhd = RectNbhdMaker(c(topl_row = 2, topl_col = 2, botr_row = 5, botr_col = 5)),
                   agents = list(
                        standard = list(
                             n = 100,
                             area_pmf = numeric()
                             ),
                        attack = list(
                             n = 0,
                             area_pmf = numeric()
                             )
                        )
                   ),
     hgrp2 = list( nbhd = RectNbhdMaker(c(topl_row = 4, topl_col = 6, botr_row = 6, botr_col = 8)),
                   agents = list(
                        standard = list(
                             n = 100,
                             area_pmf = numeric()
                             ),
                        attack = list(
                             n = 0,
                             area_pmf = numeric()
                             )
                        )
                   ),
     hgrp3 = list( nbhd = RectNbhdMaker(c(topl_row = 8, topl_col = 4, botr_row = 11, botr_col = 7)),
                   agents = list(
                        standard = list(
                             n = 100,
                             area_pmf = numeric()
                        ),
                        attack = list(
                             n = 0,
                             area_pmf = numeric()
                        )
                   )
     )    
)

HOSTILES_NAMES <- names(HOSTILES)

area <- character()
row <- numeric()
col <- numeric()
for( hgn in HOSTILES_NAMES ){
     nbhd <- HOSTILES[[hgn]]$nbhd
     area <- c(area, rep(hgn, times = nrow(nbhd)))
     row <- c(row, nbhd[,'row'])
     col <- c(col, nbhd[,'col'])
}
GRID_LAYOUT <- data.frame( area = area, row = row, col = col)
rm( nbhd, area, row, col ) #cleanup

ether <- GRID
for( hgn in HOSTILES_NAMES ){
     ether <- ether[apply( ether, MARGIN = 1, FUN = function(x){ !all(is.element(x,GRID_LAYOUT[[hgn]])) } ),]
}
if(nrow(ether) > 0){
     GRID_LAYOUT <- rbind(GRID_LAYOUT, 
                          data.frame( area = rep('ether', times = nrow(ether)), row = ether[,'row'], col = ether[,'col'] ))
}
rm( ether ) #cleanup
AREA_NAMES <- unique(as.character(GRID_LAYOUT$area))
AREA_NAMES

group <- character()
agent_class <- character()
for( hgn in HOSTILES_NAMES ){
     agent_types <- names( HOSTILES[[hgn]]$agents )
     for( atype in agent_types ){
          n <- HOSTILES[[hgn]]$agents[[atype]]$n
          group <- c(group, as.character(rep(hgn, times = n)))
          agent_class <- c(agent_class, as.character(rep(atype, times = n)))
     }
}
HOSTILES_DF <- data.frame( group = group, agent_class = agent_class)

apmf <- AreaPMFInit( group_names = HOSTILES_NAMES, 
                     area_names = AREA_NAMES,
                     grid_layout = GRID_LAYOUT,
                     within_grp_p = 0.90,
                     opposing_grp_p = 0.01,
                     ether_p = 0.09 )
HOSTILES_DF <- merge( HOSTILES_DF, as.data.frame( apmf ), by.x = "group", by.y = 0 )
v <- rviolence( n = nrow(HOSTILES_DF) )
HOSTILES_DF$violence_fac <- v
rm( group, agent_class, agent_types, hgn, atype, n, apmf, v ) #cleanup

#Look at 15 random agents from HOSTILES_DF
HOSTILES_DF[sort(sample(1:nrow(HOSTILES_DF), size = 15)), ]

AUTHORITY_DF <- data.frame()

ANGER_MATRIX <- AngerMatrixInit(names = HOSTILES_NAMES, min_anger = MIN_ANGER)



NEARBY_LIST <- ConstructNearbyList( GRID_LAYOUT, NEARBY_RADIUS )





HOSTILITY_RECORD_DF <- data.frame( time = numeric(), 
                                   perp_grp = character(), 
                                   vict_grp = character(), 
                                   row = numeric(),
                                   col = numeric() )

tot_steps <- 150

for( timestep in seq_len( tot_steps ) ){
     HOSTILES_DF <- MoveUpdater( HOSTILES_DF, area_names = AREA_NAMES, grid_layout = GRID_LAYOUT)
     print( ggplot( data = HOSTILES_DF[sample(rownames(HOSTILES_DF), size = nrow(HOSTILES_DF)), ], 
                    aes( x = col, y = -row, col=group ),   ) +
                 geom_point(alpha = 0.33) +
                 coord_cartesian( xlim = c(0,13), ylim=c(-13,0)) 
     )
     HOSTILITY_RECORD_DF <- HostilityUpdater(  current_timestep = timestep,
                                               hostility_record_df = HOSTILITY_RECORD_DF,
                                               hostiles_df = HOSTILES_DF, 
                                               authority_df = matrix(),
                                               nearby_list = NEARBY_LIST,
                                               anger_matrix = ANGER_MATRIX,
                                               lam=0.002)
     HOSTILITY_MATRIX <- HostilityMatrixMaker(group_names = HOSTILES_NAMES, 
                                              hostility_record_df = HOSTILITY_RECORD_DF, 
                                              timestep = timestep)
     

     ANGER_MATRIX <- AngerMatrixUpdater(hostility_matrix = HOSTILITY_MATRIX, 
                                        anger_matrix = ANGER_MATRIX, 
                                        min_anger = MIN_ANGER )
     print( sprintf("Time step: %d",timestep) )
     print("Hostility Matrix:")
     print(HOSTILITY_MATRIX)
     print("Anger Matrix:")
     print(ANGER_MATRIX)
     print("--------------------------------")
     
}

shots <- HOSTILITY_RECORD_DF %>% group_by( time ) %>% summarise( n = n() )
tdf <- data.frame(time = 1:tot_steps, n = 0)
tdf[shots$time,'n'] <- shots$n
shots <- tdf
ggplot( shots, aes(x = time, y = n) ) + geom_line()