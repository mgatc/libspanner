#ifndef LIBSPANNER_DN97_H
#define LIBSPANNER_DN97_H

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <vector>

#include "libspanner/greedy/types.h"

namespace spanner {

    namespace dn97 {

        void SINGLE_SOURCE(const PointContainer& V,
                const AdjacencyListDense& G,   // input graph
                           index_t i,                       // source node
                           double R1,                    // small radius
                           double R2,                    // large radius (R1 <= R2)
                           std::list<index_t>& C1,               // nodes within dist <= R1
                           std::list<index_t>& C2,               // nodes within R1 < dist <= R2
                           std::list<index_t>& C3,               // nodes having updated labels
                           std::vector<double>& ShortestPaths,     // distance from source
                           std::vector<index_t>& Parents,     // predecessor node
                           bool addToC1 = true
        )
        {
//            index_t v;
//            Edge e;
//
//            dist[s] = 0.0;
//            PQ.emplace(0.0, s);
//            C3.push_back(s);
//
//            while (!PQ.empty())
//            {
//                index_t u = (PQ.begin())->second; // add u to S
//                double du = dist[u];
//                if (du <= R1)
//                {
//                    if (u != s) C1.push_back(u);  // add u to set C1
//                }
//                else
//               {
//                    if (du <= R2)
//                        C2.push_back(u);   // add u to set C2
//                    else
//                        break;        // large cluster radius exceeded - exit
//                }
//
//                forall_adj_edges(e, u)
//                {
//                    v = G.opposite(u, e); // makes it work for ugraphs
//                    double c = du + G[e];
//                    if (pred[v] == nil && v != s )
//                    {
//                        PQ.insert(v, c); // first message to v
//                        C3.push(v);
//                    }
//                    else if (c < dist[v])
//                        PQ.decrease_p(v, c); // better path
//                    else continue;
//                    dist[v] = c;
//                    pred[v] = e;
//                }
//            }



            typedef std::map<number_t, index_t>
                    Heap;
            typedef Heap::iterator
                    HeapHandle;

            const index_t n = V.size();
            auto startPoint = V.at(i);

            Heap open;
            std::unordered_map<index_t, HeapHandle> handleToHeap(n);
            handleToHeap[i] = open.emplace(0, i).first;

            ShortestPaths[i] = 0;
            std::unordered_set<index_t> C3Helper;
            C3Helper.insert(i);

            auto current = open.begin(); // initialize current vertex to start
            index_t u_index = current->second;
            auto currentPoint = startPoint;
            auto neighborPoint = currentPoint;
            number_t newScore = 0;

            do {
                current = open.begin();

                u_index = current->second;
                currentPoint = V[u_index];

                double du = ShortestPaths[u_index];

                if (addToC1 && u_index != i && du <= R1)
                    C1.push_back(u_index);  // add u to set C1
                else if (du <= R2)
                    C2.push_back(u_index);   // add u to set C2
                else
                    break;        // large cluster radius exceeded - exit




                // loop through neighbors of current
                for (index_t neighbor : G.at(u_index)) {
                    neighborPoint = V[neighbor];
                    C3Helper.insert(neighbor);

                    newScore = du + getDistance(currentPoint, neighborPoint);

                    if (newScore < ShortestPaths[neighbor]) {
                        Parents[neighbor] = u_index;
                        ShortestPaths[neighbor] = newScore;

                        if (contains(handleToHeap, neighbor)) {
                            open.erase(handleToHeap[neighbor]);
                        }
                        handleToHeap[neighbor] = open.emplace(ShortestPaths[neighbor], neighbor).first;
                    }
                }

                open.erase(current);
            } while (!open.empty());

            std::copy(C3Helper.begin(), C3Helper.end(), std::back_inserter(C3));

        } // end SINGLE_SOURCE






    }


    void ClusterGraph(const greedy::input_t& P, greedy::output_t& out,
                      double W = 1.0, double delta = 0.8) {
        index_t M = 5; // ? What is this? Num edges. Num edges of what? Complete graph of P?
        index_t D = 2; // dimension
        index_t N = P.size();

        std::vector<Edge> L(M);

        double sp_length = 0.0;

        AdjacencyListDense G(N);


        for(Edge e : L) {
            G[e.first].insert(e.second);
            G[e.second].insert(e.first);
        }

        /**********************************************
         Construct the cluster graph CG for specific radius r
        ***********************************************/

        AdjacencyListDense CG(N); // cluster graph

        double d;


        std::vector<index_t> LCCenters;
        std::vector<index_t> LNC(N);
        std::vector<Edge> LInC(N*N);
        std::vector<Edge> LOutC1(N*N);
        std::vector<Edge> LOutC2(N*N);
        int InM, OutM1, OutM2;

        for (int i = 0; i < N; i++)
        {
//                v = CG.new_node(i);
            LNC[i] = i;
        }
        LCCenters.clear();

        InM = 0; OutM1 = 0; OutM2 = 0;
        std::vector<bool>   Mark(N, false);
        std::vector<bool>   IsClusterCenter(N, false);
        std::vector<double> dist(N, INF);
        std::vector<index_t>   pred(N);
        std::vector< std::vector< std::pair<index_t, double> > > cluster_neighbours(N);
        std::list<index_t> C1;
        std::list<index_t> C2;
        std::list<index_t> C3;

        for(index_t v=0; v<N; ++v) {
            if(!Mark[v]) {
                LCCenters.push_back(v);
                IsClusterCenter[v] = true;

                dn97::SINGLE_SOURCE(P, G, v, delta * W, W, C1, C2, C3, dist, pred);
                index_t u;
                // add intra-cluster edges
                while (!C1.empty())
                {
                    u = C1.back(); // changed from pop to pop_back, maybe an issue down the road
                    C1.pop_back();
                    Edge e = std::make_pair(LNC[u], LNC[v]);
                    CG[e.first].insert(e.second);
                    CG[e.second].insert(e.first);
                    LInC[InM++] = e;
                    Mark[u] = true;
                }

                // save information about inter-cluster neighbours
                while (!C2.empty())
                {
                    u = C2.back();
                    C2.pop_back();
                    std::pair<index_t, double> cn(v, dist[u]);
                    cluster_neighbours[u].push_back( cn );
                }

                // reset node labels
                while (!C3.empty())
                {
                    u = C3.back();
                    C3.pop_back();
                    pred[u] = SIZE_T_MAX;
                    dist[u] = INF;
                }
            }
        }

        // add inter-cluster edges of type I
        // i.e., connect cluster centers within distance W from each other

        for(index_t v : LCCenters)
        {
//                std::pair<index_t, double> cn;
            for( auto cn : cluster_neighbours[v])
            {
                index_t u = cn.first;   // node
                double d = cn.second;  // shortest path distance
                assert (IsClusterCenter[u]);
                assert (IsClusterCenter[v]);
                Edge e = std::make_pair(LNC[u], LNC[v]);
                CG[e.first].insert(e.second);
                CG[e.second].insert(e.first);
                LOutC1[OutM1++] = e;
            }
        }

        // add inter-cluster edges of type II
        // i.e., for each edge in G, add an edge between the respective cluster centers
        forall_edges(e, G)
        {
            node cu, cv;
            Edge eu, ev;
            list<Edge> clu, clv;

            u = LNC[G[G.source(e)]]; // first endpoint in GC
            v = LNC[G[G.target(e)]]; // second endpoint in GC

            // get list of cluster centers that node u belongs to
            if (IsClusterCenter[G.source(e)])
                clu.append( (Edge) 0); // dummy edge
            else
                forall_out_edges(ee, u) clu.append(ee);

            // get list of cluster centers that node v belongs to
            if (IsClusterCenter[G.target(e)])
                clv.append( (Edge) 0); // dummy edge
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
    }

    void DN97(const greedy::input_t& in, greedy::output_t& out,
              double t = 1.5, double delta = 0.2, double alpha = 2.0) {
        using namespace dn97;

        const auto& P = in;
        index_t N = P.size();
        index_t M;
        index_t D = 2;

        std::vector<Edge> clusterGraphEdges;
        ClusterGraph(P, clusterGraphEdges);


        AdjacencyListDense G(N);
        std::transform(clusterGraphEdges.begin(), clusterGraphEdges.end(), std::inserter(G,G.end()),
            [&G](const Edge& e) {
                G[e.first].insert(e.second);
                G[e.second].insert(e.first);
            });

        //**********************************************
        // Construct the spanner graph SG
        // *********************************************

        AdjacencyListDense SG;




        std::vector<index_t> LNS(N);

        SG.clear();
        for (int i = 0; i < N; i++)
        {
            LNS[i] = i;
        }

        // sort the edges of G by length
        std::map<Edge,double> Cost;
        std::transform(clusterGraphEdges.begin(),clusterGraphEdges.end(), std::inserter(Cost,Cost.end()),
            [&P](const Edge& e) {
                return std::make_pair(e, getDistance(P[e.first],P[e.second]));
            });

        std::vector<Edge> clusterGraphEdgesSorted(clusterGraphEdges);
        std::sort(clusterGraphEdgesSorted.begin(), clusterGraphEdgesSorted.end(),
            [&Cost] (const Edge& a, const Edge& b) {
                return Cost[a] < Cost[b];
            });

        double LMax = Cost[clusterGraphEdgesSorted.back()];
        int bucket  = 0;
        double R0   = LMax / (double) N, R = R0;

        std::vector<double> dist(N, INF);

        std::vector< std::list<Edge> > AdjEdges(N);
        std::map<Edge,bool> InSpanner;

        auto ee = clusterGraphEdgesSorted.begin();
        for (ee=ee; ee!=clusterGraphEdgesSorted.end() && (Cost[*ee] <= R0); ++ee)
        {
            SG[LNS[(*ee).first]].insert( LNS[(*ee).second]);
            SG[LNS[(*ee).second]].insert(LNS[(*ee).first]);
            InSpanner[*ee] = true;
        }
        for (ee=ee; ee!=clusterGraphEdgesSorted.end(); ++ee)
            {
                while (Cost[*ee] > R)
                    R = alpha*R;

                // cerr << "Processing new bucket .." << endl;
                // clear bucket
                for( auto v: AdjEdges)
                    v.clear();

                // collect edges from this bucket
//                while ((*ee) && (Cost[*ee] <= R))
                for (ee=ee; ee!=clusterGraphEdgesSorted.end() && (Cost[*ee] <= R); ++ee)
                {
                    // cerr << "Added to AdjEdge of " << G[G.source(ee)] << " and "
                    //      << G[G.target(ee)] << endl;
                    index_t su = LNS[ee->first];
                    index_t sv = LNS[ee->second];
                    AdjEdges[su].push_back(*ee);
                    AdjEdges[sv].push_back(*ee);
                }

                index_t u, uu, v, w;

                std::vector<bool> Mark(N,false);
                std::vector<double> dist(N,INF);
                std::list<index_t> C1;
                std::list<index_t> C2;
                std::list<index_t> C3;
                node_pq<double> PQ(SG);

                for(auto v : P ){}
                if (!Mark[v])
                {
                    C1.clear();
                    C2.clear();
                    while (!C3.empty())
                    {
                        uu = C3.back();
                        C3.pop_back();
                        dist[uu] = INF;
                    }
                    SINGLE_SOURCE(P, SG, v, delta*R, alpha*t*R, C1, C2, C3, dist, {}, true);

                    // process edges for this cluster center v
                    while (!C1.empty())
                    {
                        u = C1.back();
                        C1.pop_back();
                        Mark[u] = true;

                        for(auto e : AdjEdges[u])
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
}

#endif //LIBSPANNER_DN97_H