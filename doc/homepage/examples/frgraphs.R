
frgraphs <- function() {
  graphs <- list()
  graphs[[1]] <- graph.full(5)
  graphs[[2]] <- graph.full(4)
  graphs[[3]] <- graph.full(6)
  graphs[[4]] <- graph.full(6)
  graphs[[5]] <- graph.full(8)
  graphs[[6]] <- graph.full(2)
  graphs[[7]] <- graph.formula( A-B-C-D-A, A-F, B-G, C-H, D-I )
  graphs[[8]] <- graph.formula( A-B-C-A, D-E-F-D, A-D )
  graphs[[9]] <- graph.ring(5)
  graphs[[10]] <- graph.formula( A-B-C-A, B-D-E-B, C-E-F-C,
                                D-G-H-D, E-H-I-E, F-I-J-F )
  graphs[[11]] <- graph.formula( A-B-C-D-A, A-C, B-D,
                                A-E-F-G-A, A-F, E-G,
                                B-H-I-J-B, B-I, H-J,
                                C-K-L-M-C, C-L, K-M,
                                D-N-O-P-D, D-O, N-P )
  graphs[[12]] <- graph.formula( A-B-C-D-F-E-A, A-D, E-C, B-F)
  graphs[[13]] <- graph.formula( A-B:D:H, B-F:G, C-D:E:I, D-F, E-F:L, F-G:I:J:M,
                                G-N, H-K:L, I-K, J-M:N, K-M, L-N )
  graphs[[14]] <- simplify(graph.formula( A-B:H:L, B-A:C:K, C-B:D:M, D-C:E:L:M:N,
                                         E-D:F:G, F-E:I:L, G-E:H:N, H-A:G:J,
                                         I-F:J:K, J-H:I:M, K-B:I:N, L-A:D:F,
                                         M-C:D:J, N-D:G:K ))
  graphs[[15]] <- graph.formula( A-B-C-D-E-F-G-H-I-A,
                                J-E:A:F, K-A:F:B, L-F:B:G, M-B:G:C, N-G:C:H,
                                O-C:H:D, P-H:D:I, Q-D:I:E, R-I:E:A )
  graphs[[16]] <- graph.formula( 0-1:4:8:11:15:19:22, 1-0:2:6:24, 2-1:3:5:23:29:30,
                                3-2, 4-0:5, 5-2:4:6:9, 6-1:5:8:31:7:30, 7-6:30:31,
                                8-0:6:9:13, 9-5:8:12:32:10:31, 10-9, 11-0:12,
                                12-9:11:26:13, 13-8:12:27:33:14:32, 14-13:32:33,
                                15-0:16, 16-15:20:17:26, 17-16:19:35:18:34:27,
                                18-17:34:35, 19-0:24:20:17, 20-19:23:36:21:35:16,
                                21-20, 22-0:23, 23-22:2:24:20, 24-23:1:29:25:36:19,
                                25-24:29:36, 26-16:12:33:27:34, 27-26:13:28:17,
                                28-27, 29-2:24:25, 30-2:6:7, 31-6:7:9,
                                32-9:13:14, 33-13:14:26, 34-26:17:18, 35-17:18:20,
                                36-20:24:25 )
  graphs[[16]] <- simplify(graphs[[16]])
  graphs
}
