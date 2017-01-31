## Create Coordinate List for Rectangular Nbh inputing top-left (row,col) and bottom-right (row,col)
RectNbhdMaker <- function( tl_br ){ #tl_br = c( tlr, tlc, brr, brc )
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
     forget_fac <- 0.9
     prev_fac <- prev_anger * forget_fac
     
     hostility_fac <- 0.0
     if( n_hostilities == 1 ){
          hostility_fac <- 0.5
     }
     else if( n_hostilities == 2){
          hostility_fac <- 0.75
     }
     else if( n_hostilities >= 3){
          hostility_fac <- 1.0
     }
     
     return( max( hostility_fac + (1 - hostility_fac) * prev_fac, min_anger ) )
     
}

AngerMatrixUpdater <- function( hostility_matrix, anger_matrix ){
        n <- nrow(hostility_matrix)
        for( i in 1:n ){
                for( j in 1:n ){
                     if( i!=j ){
                          # note j,i vs i,j below: that way host&tens matrices are aligned the same: row=to, col=from
                          anger_matrix[i,j] <- AngerCalculator( hostility_matrix[j,i], anger_matrix[i,j] )
                     }
                }
        }
        return( anger_matrix )
}


TensionCalculator <- function( grp1, grp2, anger_matrix ){
        # given the asymmetric anger_matrix, return the symmetric tension between groups
        return( max(anger_matrix[grp1,grp2], anger_matrix[grp2,grp1]) )
}

ShootDecision <- function( perp_grp, perp_viol, vict_grp, anger_matrix, lam = 0.5 ){
     grp_tens <- TensionCalculator( perp_grp, vict_grp, anger_matrix )
     p_shoot <- 1 - exp(- lam * (perp_viol * grp_tens))
     return( sample( c(TRUE,FALSE), size = 1, prob = c(p_shoot,1-p_shoot) ) )
}

#HostilityMatrixUpdater <- function( agents_df, anger_matrix, nbr_radius ){
     # initialize new matrix to return
#     updated_hostility_matrix <- HostilityMatrixInit( hostile_group_names )
     # a function to agregate over neighbor agents
#     AgentNeighbors <- function ( agent1, agent2, radius ){
#          b1 = c(agent1['row'], agent1['col'])
#          b2 = c(agent2['row'], agent2['col'])
#          return( AreNeighbors(b1,b2,radius) )
#     }
     # no one shoots near an authority
     #p_vs_h <- agents_df$group == 'police'
     #police <- agents_df[ p_vs_h, ]
     #hostiles <- agents_df[ !p_vs_h, ]
     # this loop will remove all hostile agents from shooting if they are near police
     #not_near_police <- logical()
     #for( i in 1:nrow(police) ){
     #     p <- police[i,]
     #     for( j in 1:nrow(hostiles) ){
     #          h <- hostiles[j,]
     #          not_near_police <- c(not_near_police, !AgentNeighbors(p,h,nbr_radius))
     #     }
     #}
     #hostiles <- hostiles[ not_near_police, ]
     #at this point, every row left in hostiles is not a neighbor of police
     
#     rem_host_names <- unique( as.character(hostiles$group) )
#     for( perp_grp in rem_host_names ){
#          in_grp <- hostiles$group == perp_grp
#          perps <- hostiles[ in_grp, ]
#          victs <- hostiles[ !in_grp, ]
#          for( i in 1:nrow(perps) ){
#               perp <- perps[i, ]
               #which opposing groups are a neighbor of this perp
#               nb <- logical()
#               for( j in 1:nrow(victs) ){
#                    vict <- victs[j, ]
#                    nb <- c(nb, AgentNeighbors(perp,vict,nbr_radius))
#               }
#               nbhs <- unique( as.character( victs[ nb, ]$group ) )
#               if( length(nbhs) != 0) {
#                    for( vict_grp in nbhs ){
#                         if( ShootDecision(perp_grp = perp_grp, 
#                                           perp_viol = perp$violence_fac, 
#                                           vict_grp = vict_grp, 
#                                           anger_matrix = anger_matrix) ){
#                              updated_hostility_matrix[perp$group,vict_grp] = update_hostility_matrix[perp_grp,vict_grp] + 1
#                         }
#                    }
#               }
#          }
#     }
#     return( updated_hostility_matrix )
#}

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
                     within_grp_p = 0.8,
                     opposing_grp_p = 0.05,
                     ether_p = 0.15 )
HOSTILES_DF <- merge( HOSTILES_DF, as.data.frame( apmf ), by.x = "group", by.y = 0 )
v <- rviolence( n = nrow(HOSTILES_DF) )
HOSTILES_DF$violence_fac <- v
rm( group, agent_class, agent_types, hgn, atype, n, apmf, v ) #cleanup

#Look at 15 random agents from HOSTILES_DF
HOSTILES_DF[sort(sample(1:nrow(HOSTILES_DF), size = 15)), ]

require(ggplot2)
for( i in 1:50 ){
     HOSTILES_DF <- MoveUpdater( HOSTILES_DF, area_names = AREA_NAMES, grid_layout = GRID_LAYOUT)
     print(ggplot( data = HOSTILES_DF, aes( x = col, y = -row, col=group )  ) + geom_point(alpha = 0.33))
     Sys.sleep( time = 1 )
}

