#ifndef LIBSPANNER_DN97_H
#define LIBSPANNER_DN97_H

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "libspanner/greedy/types.h"

namespace spanner {

    namespace dn97 {

//                           Metric definition
// ===================================================================
// The Euclidean (L_2) metric is used.
// ====================================================================

#define SWAP(a,b) { help = *a; *a = *b; *b = help; }
#define MAX(a,b)  ((a > b) ? (a) : (b))
#define MIN(a,b)  ((a < b) ? (a) : (b))

        struct simple_edge
        {
            int    idx1;
            int    idx2;
            double dist;
        };

        struct sp_statistics {
            long    margin_dist_count;
            long    point_dist_count;

            float   fair_split_sort_cpu;
            float   fair_split_tree_cpu;

            float   wspd_cpu;

            long    number_of_pairs;
            long    overlapping_pairs;
            float   sparse_graph_cpu;

            long    closest_total;
            long    closest_recursive1;
            long    closest_recursive2;

            long    pq_insertions;
            long    pq_findmins;
            long    pq_max_size;

            long    kruskal_edges1;
            long    kruskal_edges2;
            float   kruskal_select_cpu;
            float   kruskal_sort_cpu;
            float   kruskal_cpu;

            double  tree_length;
            double  sp_length;
            float   total_cpu;
        };

        extern struct sp_statistics SpStats;

// ----------------------------------------------------------------------
// Compute L_p distance between two points
// References into P array
// ----------------------------------------------------------------------
        static double lp_distance(double* P, int ndim, int idx1, int idx2)
        {
            double dist = 0.0, d;
            int dim;

            // Euclidean distance between two points
            for (dim = 0; dim < ndim; dim++)
            {
                d = P[idx1++] - P[idx2++];
                dist += d*d;
            }
            dist = sqrt(dist);

            return dist;
        }

        template<class NT>
        class cmp_edgesA : public leda_cmp_base<edge> {
            const edge_array<NT>* edge_cost;
        public:
            cmp_edgesA(const edge_array<NT>& cost) : edge_cost(&cost) {}

            int operator()(const edge& x, const edge& y) const
            { return compare((*edge_cost)[x],(*edge_cost)[y]); }
        };


        static char *usage =
                "\nUsage: clusterspanner\n"
                "Compute a cluster spanner and show on the screen.\n"
                "The graph is read from standard input (e.g., from another program).\n";

        struct sp_statistics SpStats;
        bool DEBUG = false;

        void SINGLE_SOURCE(const GRAPH<int,double>& G,   // input graph
                           node s,                       // source node
                           double R1,                    // small radius
                           double R2,                    // large radius (R1 <= R2)
                           list<node>& C1,               // nodes within dist <= R1
                           list<node>& C2,               // nodes within R1 < dist <= R2
                           list<node>& C3,               // nodes having updated labels
                           node_array<double>& dist,     // distance from source
                           node_array<edge>&   pred,     // predecessor node
                           node_pq<double>& PQ           // priority queue
        )
        {
            node v;
            edge e;

            dist[s] = 0;
            PQ.insert(s, 0);
            C3.push(s);

            while (!PQ.empty())
            {
                node u = PQ.del_min(); // add u to S
                double du = dist[u];
                if (du <= R1)
                {
                    if (u != s) C1.push(u);  // add u to set C1
                }
                else
                {
                    if (du <= R2)
                        C2.push(u);   // add u to set C2
                    else
                        break;        // large cluster radius exceeded - exit
                }

                forall_adj_edges(e, u)
                {
                    v = G.opposite(u, e); // makes it work for ugraphs
                    double c = du + G[e];
                    if (pred[v] == nil && v != s )
                    {
                        PQ.insert(v, c); // first message to v
                        C3.push(v);
                    }
                    else if (c < dist[v])
                        PQ.decrease_p(v, c); // better path
                    else continue;
                    dist[v] = c;
                    pred[v] = e;
                }
            }

            // if (!PQ.empty()) PQ.clear(); // clear PQ - GIVES AN ERROR WHEN DELETING PQ
            while (!PQ.empty()) PQ.del_min(); // clear PQ

        } // end SINGLE_SOURCE


        void CLUSTER_GRAPH(const GRAPH<int,double>& G,
                           double delta,
                           double W,
                           GRAPH<int,double>& CG,
                           list<node> & LCCenters,
                           array<node> & LNC,
                           array<edge>& LInC,
                           array<edge>& LOutC1,
                           array<edge>& LOutC2,
                           int & InM,
                           int & OutM1,
                           int & OutM2)
        {
            // **********************************************
            // Construct the graph CG
            // *********************************************

            double d;
            node u, v, w;
            edge e, ee;

            CG.clear();
            int N = G.number_of_nodes();
            for (int i = 0; i < N; i++)
            {
                v = CG.new_node(i);
                LNC[i] = v;
            }
            LCCenters.clear();

            InM = 0; OutM1 = 0; OutM2 = 0;
            node_array<bool>   Mark(G);
            node_array<bool>   IsClusterCenter(G);
            node_array<double> dist(G);
            node_array<edge>   pred(G);
            node_array< list< two_tuple<node, double> > > cluster_neighbours(G);
            list<node> C1;
            list<node> C2;
            list<node> C3;
            node_pq<double> PQ(G);

            forall_nodes(v, G) Mark[v] = false;
            forall_nodes(v, G) IsClusterCenter[v] = false;
            forall_nodes(v, G) pred[v] = nil;
            forall_nodes(v, G) dist[v] = MAXDOUBLE;

            forall_nodes(v, G)
            {
                if (!Mark[v])
                {
                    LCCenters.append(v);
                    IsClusterCenter[v] = true;

                    SINGLE_SOURCE(G, v, delta*W, W, C1, C2, C3, dist, pred, PQ);

                    // add intra-cluster edges
                    while (!C1.empty())
                    {
                        u = C1.pop();
                        edge e = CG.new_edge(LNC[G[u]], LNC[G[v]], dist[u]);
                        LInC[InM++] = e;
                        Mark[u] = true;
                    }

                    // save information about inter-cluster neighbours
                    while (!C2.empty())
                    {
                        u = C2.pop();
                        two_tuple<node, double> cn(v, dist[u]);
                        cluster_neighbours[u].append( cn );
                    }

                    // reset node labels
                    while (!C3.empty())
                    {
                        u = C3.pop();
                        pred[u] = nil;
                        dist[u] = MAXDOUBLE;
                    }
                }
            }

            // add inter-cluster edges of type I
            // i.e., connect cluster centers within distance W from each other

            forall(v, LCCenters)
            {
                two_tuple<node, double> cn;
                forall(cn, cluster_neighbours[v])
                {
                    u = cn.first();   // node
                    d = cn.second();  // shortest path distance
                    assert (IsClusterCenter[u]);
                    assert (IsClusterCenter[v]);
                    edge e = CG.new_edge(LNC[G[v]], LNC[G[u]], d);
                    LOutC1[OutM1++] = e;
                }
            }

            // add inter-cluster edges of type II
            // i.e., for each edge in G, add an edge between the respective cluster centers
            forall_edges(e, G)
            {
                node cu, cv;
                edge eu, ev;
                list<edge> clu, clv;

                u = LNC[G[G.source(e)]]; // first endpoint in GC
                v = LNC[G[G.target(e)]]; // second endpoint in GC

                // get list of cluster centers that node u belongs to
                if (IsClusterCenter[G.source(e)])
                    clu.append( (edge) 0); // dummy edge
                else
                    forall_out_edges(ee, u) clu.append(ee);

                // get list of cluster centers that node v belongs to
                if (IsClusterCenter[G.target(e)])
                    clv.append( (edge) 0); // dummy edge
                else
                    forall_out_edges(ee, v) clv.append(ee);

                //forall(ee, clu) { if (ee != 0) cerr << CG[CG.target(ee)] << " "; } cerr << endl;
                //forall(ee, clv) { if (ee != 0) cerr << CG[CG.target(ee)] << " "; } cerr << endl;

                forall(eu, clu)
                forall(ev, clv)
                {
                    if (eu != 0) cu = CG.target(eu); else cu = u;
                    if (ev != 0) cv = CG.target(ev); else cv = v;
                    if (cu != cv) // only distinct node clusters
                    {
                        //cerr << "cluster " << CG[cu] << " and " << CG[cv] << endl;

                        // are cluster centers already connected? if so, then get edge
                        edge uve = (edge) 0;
                        forall_inout_edges(ee, cu)
                        if (cv == CG.opposite(cu, ee)) { uve = ee; break; }

                        // add edge if necessary
                        if (uve == (edge) 0)
                        {
                            uve = CG.new_edge(LNC[G[cu]], LNC[G[cv]], MAXDOUBLE);
                            LOutC2[OutM2++] = uve;
                            if (DEBUG) cerr << "edge added!!\n";
                        }

                        // update distance

                        double newd = G[e];
                        if (!IsClusterCenter[G.source(e)])
                            newd += CG[eu]; // add distance to cluster center
                        if (!IsClusterCenter[G.target(e)])
                            newd += CG[ev]; // add distance to cluster center

                        if (CG[uve] > newd) {
                            if (DEBUG) cerr << "distance updated!!\n";
                            CG[uve] = newd;
                        }
                    }
                }
            }

        } // end CLUSTER_GRAPH


        int clusterGraphMain(int argc, char *argv[])
        {
            int N, M, D, param = 1;
            double W = 1.0;
            double Delta = 0.8;		/* Delta value for cluster */

            while (param < argc)
            {
                if ((!strcmp(argv[param],"-h")) || (!strcmp(argv[param],"-H")))
                    // Usage
                { fprintf(stderr, "%s\n", usage); return 0; }
                if ((!strcmp(argv[param],"-w")) || (!strcmp(argv[param],"-W")))
                    // Cluster radius
                    W = atof(argv[++param]);
                if ((!strcmp(argv[param],"-l")) || (!strcmp(argv[param],"-L")))
                    // Delta
                    Delta = atof(argv[++param]);
                if ((!strcmp(argv[param],"-d")) || (!strcmp(argv[param],"-D")))
                    // turn on DEBUG
                    DEBUG = true;
                param++;
            }

            /**********************************************
                      Read graph G from stdin
            ***********************************************/

            cin >> D;
            cin >> N;
            cin >> M;
            double *P = new double[N*D];
            for (int i = 0; i < N; i++)
            {
                for (int dim = 0; dim < D; dim++)
                    cin >> P[i*D + dim];
            }

            simple_edge* Edges = new simple_edge[M];
            for (int i = 0; i < M; i++)
                cin >> Edges[i].idx1 >> Edges[i].idx2 >> Edges[i].dist;

            float total_cpu = used_time();
            double sp_length = 0.0;

            /**********************************************
                Construct the graph G (only works for 2d)
            ***********************************************/


            GRAPH<int, double> G; // declared globally
            G.clear();
            node v, w;
            edge e, f;

            if (D == 2)
            {
                // create the input graph G
                array<node> LN(N);
                for (int i = 0; i < N; i++)
                {
                    v = G.new_node(i);
                    LN[i] = v;
                }

                array<edge> LE(2*M);
                for (int i = 0; i < M; i++)
                {
                    e = G.new_edge(LN[Edges[i].idx1], LN[Edges[i].idx2]);
                    LE[i] = e;
                    e = G.new_edge(LN[Edges[i].idx2], LN[Edges[i].idx1]);
                    LE[i+M] = e;
                }

                edge_array<double> Cost(G);
                forall_edges(e,G)
                {
                    double temp = lp_distance(P, D, D*G[G.source(e)], D*G[G.target(e)]);
                    Cost[e] = temp; G[e] = temp;
                }

                // sort the edges of G by length
                list<edge> L;
                L.clear();
                for (int i=0; i < M; i++)
                    L.append(LE[i]);
                cmp_edgesA<double> cmp(Cost);
                L.sort(cmp);

                /**********************************************
                 Construct the cluster graph CG for specific radius r
                ***********************************************/

                GRAPH<int, double> CG; // cluster graph
                list<node> LCCenters;
                array<node> LNC(N);
                array<edge> LInC(N*N);
                array<edge> LOutC1(N*N);
                array<edge> LOutC2(N*N);
                int InM, OutM1, OutM2;

                CLUSTER_GRAPH(G, Delta, W, CG, LCCenters, LNC,
                              LInC, LOutC1, LOutC2, InM, OutM1, OutM2);

                total_cpu = used_time(total_cpu);

                // output CG for showgraph
                cout << D << endl << N << endl << InM + OutM1 + OutM2 << endl;
                cout << LCCenters.length() << endl << Delta*W << endl;
                cout << InM << endl << OutM1 << endl << OutM2 << endl;
                for (int i = 0; i < N; i++)
                    cout << P[i*D] << ' ' << P[i*D+1] << endl;
                for (int i = 0; i < InM; i++)
                    cout << CG[CG.source(LInC[i])] << ' ' << CG[CG.target(LInC[i])]
                         << ' ' << CG[LInC[i]] << endl;
                for (int i = 0; i < OutM1; i++)
                    cout << CG[CG.source(LOutC1[i])] << ' ' << CG[CG.target(LOutC1[i])]
                         << ' ' << CG[LOutC1[i]] << endl;
                for (int i = 0; i < OutM2; i++)
                    cout << CG[CG.source(LOutC2[i])] << ' ' << CG[CG.target(LOutC2[i])]
                         << ' ' << CG[LOutC2[i]] << endl;
                forall(v,LCCenters)
                cout << G[v] << endl;

                // summary of results
                fprintf(stderr,
                        "ClusterGraph:  Edges(In) \n");
                fprintf(stderr,
                        "@6              %1d\n", M);
                fprintf(stderr,
                        "ClusterGraph:  Edges(Out) IntraEdges Inter1Edges Inter2Edges \n");
                fprintf(stderr,
                        "@7              %1d        %d         %d          %d\n",
                        InM + OutM2 + (OutM1/2), InM, OutM1/2, OutM2);

                fprintf(stderr, "ClusterGraph:  NumCenters Length TotalCPU NormTotCPU\n");
                fprintf(stderr, "@8              %1d        %.6f   %.2f   %.2f\n",
                        LCCenters.length(), sp_length, total_cpu,
                        1000000.0 * total_cpu / ((double) N * log10( (double) N)));
            }
        }


        static char *usage =
                "\nUsage: clusterspanner\n"
                "Compute a cluster spanner and show on the screen.\n"
                "The graph is read from standard input (e.g., from another program).\n";

        struct sp_statistics SpStats;

        double delta = 0.2;
        double alpha = 2;
        int number_i, number_I, number_II;

        void SINGLE_SOURCE(const GRAPH<int,double>& G,     // input graph
                           node s,                         // source node
                           double R1,                      // small radius
                           double R2,                      // large radius (R1 <= R2)
                           list<node>& C1,                 // nodes within dist <= R1
                           list<node>& C2,                 // nodes within R1 < dist <= R2
                           list<node>& C3,                 // nodes having updated labels
                           node_array<double>& dist,       // distance from source
                           node_pq<double>& PQ,            // priority queue
                           bool AddtoC1                    // add to C1 or not
        )
        {
            node v;
            edge e;

            dist[s] = 0;
            PQ.insert(s, 0);
            C3.push(s);

            while (!PQ.empty())
            {
                node u = PQ.del_min(); // add u to S
                double du = dist[u];
                if (du <= R1)
                {
                    if (AddtoC1) C1.push(u);  // add u to set C1
                }
                else
                {
                    if (du <= R2)
                        C2.push(u);   // add u to set C2
                    else
                        break;        // large cluster radius exceeded - exit
                }

                forall_adj_edges(e, u)
                {
                    v = G.opposite(u, e); // makes it work for ugraphs
                    double c = du + G[e];
                    if (dist[v] == MAXDOUBLE && v != s )
                    {
                        PQ.insert(v, c); // first message to v
                        C3.push(v);
                    }
                    else if (c < dist[v])
                        PQ.decrease_p(v, c); // better path
                    else continue;
                    dist[v] = c;
                }
            }

            // if (!PQ.empty()) PQ.clear(); // clear PQ - GIVES AN ERROR WHEN DELETING PQ
            while (!PQ.empty()) PQ.del_min(); // clear PQ

        } // end SINGLE_SOURCE


        void CLUSTER_SPANNER(GRAPH<int,double>& G,      // input graph
                             double t,                  // stretch factor
                             GRAPH<int,double>& SG      // resulting spanner graph
        )
        {
            node v, w;
            edge e, ee;

            // **********************************************
            // Construct the graphs SG
            // *********************************************

            int N = G.number_of_nodes();
            array<node> LNS(N);

            SG.clear();
            for (int i = 0; i < N; i++)
            {
                v = SG.new_node(i);
                LNS[i] = v;
            }

            edge_array<double> Cost(G);
            forall_edges(e, G)
            Cost[e] = G[e];

            // sort the edges of G by length
            G.sort_edges(Cost);

            double LMax = G[G.last_edge()];
            int bucket  = 0;
            double R0   = LMax / (double) N, R = R0;

            node_array<double> dist(SG, MAXDOUBLE);
            node_pq<double> PQ(SG);

            node_array< list<edge> > AdjEdges(SG);
            edge_array<bool> InSpanner(G,false);

            ee = G.first_edge();
            while ((ee) && (G[ee] <= R0))
            {
                SG.new_edge(LNS[G[G.source(ee)]], LNS[G[G.target(ee)]], G[ee]);
                SG.new_edge(LNS[G[G.target(ee)]], LNS[G[G.source(ee)]], G[ee]);
                // cerr << "Adding short edges .." << endl;
                InSpanner[ee] = true;
                ee = G.succ_edge(ee);
            }
            while (ee)
            {
                while (G[ee] > R)
                    R = alpha*R;

                // cerr << "Processing new bucket .." << endl;
                // clear bucket
                forall_nodes(v,SG)
                AdjEdges[v].clear();

                // collect edges from this bucket
                while ((ee) && (G[ee] <= R))
                {
                    // cerr << "Added to AdjEdge of " << G[G.source(ee)] << " and "
                    //      << G[G.target(ee)] << endl;
                    node su = LNS[G[G.source(ee)]];
                    node sv = LNS[G[G.target(ee)]];
                    AdjEdges[su].append(ee);
                    AdjEdges[sv].append(ee);
                    ee = G.succ_edge(ee);
                }

                node u, uu, v, w;

                node_array<bool>   Mark(SG,false);
                node_array<double> dist(SG,MAXDOUBLE);
                list<node> C1;
                list<node> C2;
                list<node> C3;
                node_pq<double> PQ(SG);

                forall_nodes(v, SG)
                if (!Mark[v])
                {
                    if (!C1.empty()) C1.clear();
                    if (!C2.empty()) C2.clear();
                    while (!C3.empty())
                    {
                        uu = C3.pop();
                        dist[uu] = MAXDOUBLE;
                    }
                    SINGLE_SOURCE(SG, v, delta*R, alpha*t*R, C1, C2, C3, dist, PQ, true);

                    // process edges for this cluster center v
                    while (!C1.empty())
                    {
                        u = C1.pop();
                        Mark[u] = true;

                        forall(e,AdjEdges[u])
                        {
                            if (InSpanner[e]) continue;
                            if (SG[u] == G[G.source(e)])
                                w = LNS[G[G.target(e)]];
                            else w = LNS[G[G.source(e)]];
                            // cerr << "Checking edge " << SG[u] << " to " << SG[w] << endl;

                            if ((dist[w] == MAXDOUBLE) || (dist[u] + dist[w] > t * G[e]))
                            {
                                // cerr << "Adding edge " << SG[u] << " to " << SG[w] << endl;
                                // add to SG
                                SG.new_edge(u, w, G[e]);
                                SG.new_edge(w, u, G[e]);
                                InSpanner[e] = true;

                                // redo SINGLE_SOURCE after clearing data structures
                                C2.clear();
                                while (!C3.empty())
                                {
                                    uu = C3.pop();
                                    dist[uu] = MAXDOUBLE;
                                }
                                SINGLE_SOURCE(SG, v, delta*R, alpha*t*R,
                                              C1, C2, C3, dist, PQ, false);

                            } // if
                        } // forall(e)
                    } // while
                } // if (!Mark)
            } // while (e)
        }


        int clusterSpannerMain(int argc, char *argv[])
        {
            int N, M, D, param = 1;
            double t = 1.5;

            while (param < argc)
            {
                if ((!strcmp(argv[param],"-h")) || (!strcmp(argv[param],"-H")))  // Usage
                { fprintf(stderr, "%s\n", usage); return 0; }

                if (!strcmp(argv[param],"-t"))  // Stretch factor
                    t = atof(argv[++param]);

                if (!strcmp(argv[param],"-d"))  // Delta
                    delta = atof(argv[++param]);

                if (!strcmp(argv[param],"-a"))  // Alpha
                    alpha = atof(argv[++param]);

                param++;
            }

            /**********************************************
                      Read graph G from stdin
            ***********************************************/

            cin >> D;
            cin >> N;
            cin >> M;
            double *P = new double[N*D];
            for (int i = 0; i < N; i++)
            {
                for (int dim = 0; dim < D; dim++)
                    cin >> P[i*D + dim];
            }

            simple_edge* Edges = new simple_edge[M];
            for (int i = 0; i < M; i++)
                cin >> Edges[i].idx1 >> Edges[i].idx2 >> Edges[i].dist;

            float total_cpu = used_time();
            double sp_length = 0.0;

            /**********************************************
                Construct the graph G (only works for 2d)
            ***********************************************/

            GRAPH<int, double> G;
            G.clear();
            node v, w;
            edge e, f;

            if (D != 2)
            {
                cerr << "clusterspanner: Cannot deal with dimension higher than 2" << endl;
                exit(0);
            }

            // create the input graph G
            array<node> LN(N);             // map index i -> vertex of G
            for (int i = 0; i < N; i++)
            {
                v = G.new_node(i);
                LN[i] = v;
            }

            for (int i = 0; i < M; i++)
            {
                G.new_edge(LN[Edges[i].idx1], LN[Edges[i].idx2]);
                G.new_edge(LN[Edges[i].idx2], LN[Edges[i].idx1]);
            }

            forall_edges(e, G)
            G[e] = lp_distance(P, D, D*G[G.source(e)], D*G[G.target(e)]);

            //**********************************************
            // Construct the spanner graph SG
            // *********************************************

            GRAPH<int, double> SG;
            CLUSTER_SPANNER(G, t, SG);

            total_cpu = used_time(total_cpu);

            // output CG for showgraph
            cout << D << endl << N << endl << SG.number_of_edges() << endl;
            for (int i = 0; i < N; i++)
                cout << P[i*D] << ' ' << P[i*D+1] << endl;
            forall_edges(e, SG)
            cout << SG[SG.source(e)] << ' ' << SG[SG.target(e)]
                 << ' ' << SG[e] << endl;

            // summary of results
            fprintf(stderr, "ClusterSpanner:  SF(reqd) NumEdges(In)\n");
            fprintf(stderr, "@2              %.2f     %1d    \n", t, M);

            fprintf(stderr, "ClusterSpanner:  NumEdges(Out) Length TotalCPU NormTotCPU\n");
            fprintf(stderr, "@3              %1d           %.6f   %.2f   %.2f\n",
                    SG.number_of_edges()/2, sp_length, total_cpu,
                    1000000.0 * total_cpu / ((double) N * log10( (double) N)));
        }

    }
}

#endif //LIBSPANNER_DN97_H