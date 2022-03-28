#ifndef LIBSPANNER_DN97_H
#define LIBSPANNER_DN97_H

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <vector>

#include "../greedy/types.h"

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
        std::vector<std::unordered_map<index_t,double>> CGdist(N);

        double d;


        std::vector<index_t> LCCenters;
        std::vector<index_t> LNC(N);
        std::vector<Edge> LInC(N*N);
        std::vector<Edge> LOutC1(N*N);
        std::vector<Edge> LOutC2(N*N);
        int InM, OutM1, OutM2;

        for (index_t i = 0; i < N; i++)
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
        Edge infEdge = std::make_pair(SIZE_T_MAX,SIZE_T_MAX);


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
                    double edgeLength = getEdgeLength(e,P);

                    CG[e.first].insert(e.second);
                    CGdist[e.first][e.second] = edgeLength;

                    CG[e.second].insert(e.first);
                    CGdist[e.second][e.first] = edgeLength;
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
                double edgeLength = getEdgeLength(e,P);
                CGdist[e.first][e.second] = edgeLength;
                CGdist[e.second][e.first] = edgeLength;

                CG[e.first].insert(e.second);
                CG[e.second].insert(e.first);
                LOutC1[OutM1++] = e;
            }
        }

        // add inter-cluster edges of type II
        // i.e., for each edge in G, add an edge between the respective cluster centers
//        forall_edges(e, G)
        for(auto e:L)
        {
            index_t cu, cv;
            Edge eu, ev;
            std::list<Edge> clu, clv;

            index_t u = LNC[e.first]; // first endpoint in GC
            index_t v = LNC[e.second]; // second endpoint in GC

            // get list of cluster centers that node u belongs to
            if (IsClusterCenter[e.first])
                clu.push_back(infEdge); // dummy edge
            else
//                forall_out_edges(ee, u) clu.append(ee);
                for(auto ee : G[u]) clu.emplace_back(u,ee);
            // get list of cluster centers that node v belongs to
            if (IsClusterCenter[e.second])
                clv.push_back(infEdge); // dummy edge
            else
//                forall_out_edges(ee, v) clv.append(ee);
                for(auto ee : G[v]) clv.emplace_back(v,ee);

            //forall(ee, clu) { if (ee != 0) cerr << CG[CG.target(ee)] << " "; } cerr << endl;
            //forall(ee, clv) { if (ee != 0) cerr << CG[CG.target(ee)] << " "; } cerr << endl;

//            forall(eu, clu)
            for(auto eu : clu)
            {
//                forall(ev, clv)
                for(auto ev : clv)
                {
                    if (eu.first != SIZE_T_MAX && eu.second != SIZE_T_MAX) cu = eu.second; else cu = u;
                    if (ev.first != SIZE_T_MAX && ev.second != SIZE_T_MAX) cv = ev.second; else cv = v;
                    if (cu != cv) // only distinct node clusters
                    {
                        //cerr << "cluster " << CG[cu] << " and " << CG[cv] << endl;

                        // are cluster centers already connected? if so, then get edge
                        Edge uve = std::make_pair(SIZE_T_MAX,SIZE_T_MAX);
//                        forall_inout_edges(ee, cu)
                        for(auto p:G[cu]) {
                            if (cv == p) {
                                uve = std::make_pair(cu,p);
                                break;
                            }
                        }

                        // add edge if necessary
                        if (uve == infEdge) {
                            CG[LNC[cu]].insert(LNC[cv]);
                            CG[LNC[cv]].insert(LNC[cu]);
                            double edgeLength = getEdgeLength({LNC[cu],LNC[cv]},P);
                            CGdist[LNC[cu]][LNC[cv]] = edgeLength;
                            CGdist[LNC[cv]][LNC[cu]] = edgeLength;
//                            uve = CG.new_edge(LNC[G[cu]], LNC[G[cv]], INF);
                            LOutC2[OutM2++] = uve;
//                            if (DEBUG) cerr << "edge added!!\n";
                        }

                        // update distance

                        double newd = getEdgeLength(e,P);//G[e];
                        if (!IsClusterCenter[e.first])
                            newd += getEdgeLength(eu,P); // add distance to cluster center
                        if (!IsClusterCenter[e.second])
                            newd += getEdgeLength(ev,P); // add distance to cluster center

                        if (getEdgeLength(uve,P) > newd) {
//                            if (DEBUG) cerr << "distance updated!!\n";
//                            CG[uve] = newd;
                            double edgeLength = getEdgeLength(uve,P);
                            CGdist[uve.first][uve.second] = edgeLength;
                            CGdist[uve.second][uve.first] = edgeLength;
                        }
                    }
                }
            }
        }
        for( index_t i=0; i<N; ++i) {
            std::transform(CG[i].begin(),CG[i].end(),std::back_inserter(out),
                [i](const auto& j){
                    return std::make_pair(i,j);
                });

        }
    }

    void DN97(const greedy::input_t& in, greedy::output_t& out,
              double t = 1.5, double delta = 0.2, double alpha = 2.0) {
        using namespace dn97;

        const auto& P = in;
        index_t N = P.size();

        if( N <= 1) return;

        index_t M;
        index_t D = 2;

        std::vector<Edge> clusterGraphEdges;
        ClusterGraph(P, clusterGraphEdges);

        assert(!clusterGraphEdges.empty());

        AdjacencyListDense G(N);
//        std::transform(clusterGraphEdges.begin(), clusterGraphEdges.end(), std::inserter(G,G.end()),
        for(auto e : clusterGraphEdges){
            G[e.first].insert(e.second);
            G[e.second].insert(e.first);
        }

        //**********************************************
        // Construct the spanner graph SG
        // *********************************************

        AdjacencyListDense SG;




        std::vector<index_t> LNS(N);

        SG.clear();
        for (index_t i = 0; i < N; i++)
        {
            LNS[i] = i;
        }

        // sort the edges of G by length
        std::map<Edge,double> Cost;
        std::transform(clusterGraphEdges.begin(),clusterGraphEdges.end(), std::inserter(Cost,Cost.end()),
            [&P, &Cost](const Edge& e) {
                Cost[e] = getDistance(P[e.first],P[e.second]);
                return std::make_pair(reverse_pair(e), Cost[e]);
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

                index_t u = SIZE_T_MAX,
                    uu = SIZE_T_MAX, v = SIZE_T_MAX, w = SIZE_T_MAX;

                std::vector<bool> Mark(N,false);
                std::vector<double> dist(N,INF);
                std::list<index_t> C1;
                std::list<index_t> C2;
                std::list<index_t> C3;

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
                    std::vector<index_t> tempParents(N, SIZE_T_MAX);
                    dn97::SINGLE_SOURCE(P, SG, v, delta*R, alpha*t*R, C1, C2, C3, dist, tempParents, true);

                    // process edges for this cluster center v
                    while (!C1.empty())
                    {
                        u = C1.back();
                        C1.pop_back();
                        Mark[u] = true;

                        for(auto e : AdjEdges[u])
                        {
                            if (InSpanner[e]) continue;
                            if (u == e.first)
                                w = LNS[e.second];
                            else w = LNS[e.first];
                            // cerr << "Checking edge " << SG[u] << " to " << SG[w] << endl;

                            if ((dist[w] == INF) || (dist[u] + dist[w] > t * getEdgeLength(e,P)))
                            {
                                // cerr << "Adding edge " << SG[u] << " to " << SG[w] << endl;
                                // add to SG
                                SG[u].insert(w);//, G[e]);
                                SG[w].insert(u);//.new_edge(w, u, G[e]);
                                InSpanner[e] = true;

                                // redo SINGLE_SOURCE after clearing data structures
                                C2.clear();
                                while (!C3.empty())
                                {
                                    uu = C3.back();
                                    C3.pop_back();
                                    dist[uu] = INF;
                                }
                                std::vector<index_t> tempParents2(N, SIZE_T_MAX);
                                dn97::SINGLE_SOURCE(P, SG, v, delta*R, alpha*t*R,
                                              C1, C2, C3, dist, tempParents2, false);

                            } // if
                        } // forall(e)
                    } // while
                } // if (!Mark)
            } // while (e)



    }
}

#endif //LIBSPANNER_DN97_H