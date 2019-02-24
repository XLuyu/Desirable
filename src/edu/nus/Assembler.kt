package edu.nus

import java.io.File


class Assembler(K: Int) {
    val nodes = HashMap<Pair<Long,Int>,Node>()
    val contigs = mutableListOf<String>()
    val mask = (1L shl 2*K)-1

    data class Edge(val readSet: HashSet<Long>, var seq:StringBuilder){
        // seq is store as if append after this kmer
        constructor(t: HashSet<Long>, c:Char) : this(t,StringBuilder(1).append(c))
    }
    class Node(val kmer:Long, val end:Int) {
        val edges:ArrayList<Edge> = ArrayList<Edge>(1)
        fun dest(e:Edge):Pair<Long,Int> {
            val t = (Util.decode(kmer).let{if (end==0) it.revCompl() else it}+e.seq).takeLast(Desirable.K)
            return Util.canonical(t).let { Pair(it.first,it.second xor 1) }
        }
        fun dest(i:Int):Pair<Long,Int>{
            val t = (Util.decode(kmer).let{if (end==0) it.revCompl() else it}+edges[i].seq).takeLast(Desirable.K)
            return Util.canonical(t).let { Pair(it.first,it.second xor 1) }
        }
        fun findEdgeTo(t:Node):Edge? = edges.find { dest(it)==Pair(t.kmer,t.end) }
    }
    fun Node.otherEnd() = nodes[Pair(kmer, end xor 1)]!!

    fun constructGraphFromKmerSet(kmerfile: File) {
        kmerfile.forEachLine {
            val (kmer,_) = it.split('\t')
            val ckmer = Util.canonical(kmer).first
            nodes[Pair(ckmer, 0)] = Node(ckmer, 0)
            nodes[Pair(ckmer, 1)] = Node(ckmer, 1)
        }
        println("[Kmer] ${nodes.size/2} Kmers in reads containing exclusive Kmers")
        val placeholder = HashSet<Long>()
        for ((pair,_) in nodes) if (pair.second==0) {
            val (kmer,_) = pair
            val prefix = kmer ushr 2
            val suffix = (kmer and (mask ushr 2)) shl 2
            for (c in 0..3){
                val (cp, cpr) = Util.canonical((c.toLong() shl (2*Desirable.K-2))+prefix)
                val (sc, scr) = Util.canonical(suffix+c)
                if (Pair(cp, 0) in nodes && cp<kmer) {
                    nodes[Pair(kmer, 0)]!!.edges.add(Edge(placeholder, Util.decodeTable[c xor 2]))
                    nodes[Pair(cp, 1-cpr)]!!.edges.add(Edge(placeholder, Util.decodeTable[(kmer % 4).toInt()]))
                }
                if (Pair(sc, 0) in nodes && sc<kmer) {
                    nodes[Pair(kmer,1)]!!.edges.add(Edge(placeholder, Util.decodeTable[c]))
                    nodes[Pair(sc, scr)]!!.edges.add(Edge(placeholder, Util.decodeTable[(kmer ushr (2*Desirable.K-2)).toInt() xor 2]))
                }
            }
        }
        println("[Assembly] de bruijn graph construction finished\nEdge: ${nodes.map { it.value.edges.size } .sum()}")
    }
    private fun pathContraction() {
        fun compactChain(sKmer:Long, sEnd:Int): Edge{
            var kmer = sKmer
            var end = sEnd
            val edge = Edge(HashSet<Long>(),StringBuilder())
            while (true){
                val node = nodes[Pair(kmer,end)]!!
                edge.seq.append(Util.decodeTable[(if (end==1) (kmer % 4) else kmer.ushr(2*Desirable.K-2).xor(2)).toInt()])
                if (node.edges.size!=1 || node.otherEnd().edges.size!=1) break
                val (tKmer,tEnd) = node.dest(0)
                kmer = tKmer
                end = tEnd
            }
            return edge
        }
        for ((_,node) in nodes)
            if (node.edges.size!=1 || node.otherEnd().edges.size!=1)
                for ( i in node.edges.indices) {
                    val (t,end) = node.dest(i)
                    node.edges[i] = compactChain(t,end)
                }
        for ((kmer,_) in nodes.keys.filter { it.second==0 })
            if (nodes[Pair(kmer,0)]!!.edges.size==1 && nodes[Pair(kmer,1)]!!.edges.size==1) {
                nodes.remove(Pair(kmer,0))
                nodes.remove(Pair(kmer,1))
            }
    }
    fun loadReadHopsOnGraph(file:File) {
        val k = Desirable.K
        var readID = 0L
        val seqReader = ReadFileReader(file,"Load reads hops on De Bruijn Graph")
        while (seqReader.nextRead()){
            readID += 1
            val seq = seqReader.seq!!
            val hops = (k .. seq.length).mapNotNull {
                val node = nodes[Util.canonical(seq.substring(it-k,it))]
                if (node!=null) Pair(it,node) else null
            }.toMutableList()
            if (hops.isEmpty()) continue
            if (hops.first().first!=k){
                val e = hops.first().second.edges.find { it.seq[0]== complement[seq[hops.first().first-k-1]] }
                if (e!=null) hops.add(0,Pair(-1,nodes[hops.first().second.dest(e)]!!))
            }
            if (hops.last().first!=seq.length){
                val e = hops.last().second.edges.find { it.seq[0]== seq[hops.last().first] }
                if (e!=null) hops.add(Pair(-1,nodes[hops.last().second.dest(e)]!!))
            }
            if (hops.size<=2) continue
            for (i in 1 until hops.size){
                val code = if (i==1) 3 else if (i==hops.size-1) 2 else 0
                hops[i].second.findEdgeTo(hops[i-1].second)?.readSet?.add((readID shl 2)+code)
                hops[i-1].second.otherEnd().findEdgeTo(hops[i].second.otherEnd())?.readSet?.add((readID shl 2)+(code xor 1))
            }
        }
    }

    fun getContigByTraverse() {
        val visit = HashMap<Long,Int>()
        fun dfs(sKmer:Long, sEnd:Int):StringBuilder {
            var kmer = sKmer
            var end = sEnd
            val sb = StringBuilder()
            val readSet = HashSet<Long>()
            var readCnt = 0
            while (true){
                if (visit[kmer]?:-1==readCnt) break else visit[kmer] = readCnt
                val node = nodes[Pair(kmer,end)]!!
                val edge = node.edges
                if (edge.isEmpty()) break
                val bestBranch = edge.maxBy{e-> (e.readSet.map { it ushr 2 } intersect readSet).size}!!
                readSet.addAll(bestBranch.readSet.filter { it.and(3)==2L }.map { it ushr 2 })
                val readEnd = bestBranch.readSet.filter { it.and(3)==3L }.map { it ushr 2 }
                if (bestBranch.readSet.removeIf {readSet.contains(it ushr 2)}) readCnt += 1
                readSet.removeAll(readEnd)
                sb.append(if (end==1) bestBranch.seq else UtilString.reverse(bestBranch.seq))
                val (tKmer,tEnd) = node.dest(bestBranch)
                kmer = tKmer
                end = tEnd
            }
            return sb
        }
        val leafs = arrayListOf<Pair<Pair<Long,Int>,Int>>()
        for ((pair,node) in nodes)
            if (node.edges.isEmpty()) {
                leafs.add(Pair(pair,node.otherEnd().edges.map { it.seq.length }.max()?:0 ))
            }
        leafs.sortBy { -it.second }
        for ((pair,_) in leafs)
            if (!visit.contains(pair.first)){
                val s = Util.decode(if (pair.second==1) Util.reverse(pair.first) else pair.first) + dfs(pair.first,pair.second xor 1).toString()
                contigs.add(s)
            }
        println("[Assembly] ${contigs.size} contigs generated")
    }
    fun writeContigsToFasta(file: File, lengthLimit:Int = 0){
        println("Congtig > $lengthLimit: ${contigs.count { it.length>lengthLimit }}")
        val writer = file.writer()
        for ((counter, contig) in contigs.filter { it.length>lengthLimit }.withIndex()){
            writer.write(">ExSeq${counter}_${contig.length}\n$contig\n")
        }
        writer.close()
    }
    fun run(kmer:File, reads:File, out:File,minLength:Int){
        constructGraphFromKmerSet(kmer)
        pathContraction()
        loadReadHopsOnGraph(reads)
        getContigByTraverse()
        writeContigsToFasta(out,minLength)
    }
}