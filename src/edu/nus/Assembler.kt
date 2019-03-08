package edu.nus

import java.io.File


class Assembler(var K: Int) {
    val nodes = java.util.TreeMap<Long,Node>()
    val contigs = mutableListOf<String>()
    val mask = (1L shl 2*K)-1

    data class Edge(val readSet: HashSet<Long>, var seq:StringBuilder){ // seq is store as if append right side
        constructor(t: HashSet<Long>, c:Char) : this(t,StringBuilder(1).append(c))
    }
    class Node(val kmer:Long) {
        val edges:ArrayList<Edge> = ArrayList<Edge>(1)
        fun dest(e:Edge):Long {
            val t = (Util.decode(Math.abs(kmer)).let{if (kmer>0) it.revCompl() else it}+e.seq).takeLast(Desirable.K)
            return Util.canonical(t).let { if (it.second==0) -it.first else it.first }
        }
        fun dest(i:Int) = dest(edges[i])
        fun findEdgeTo(t:Long):Edge? = edges.find { dest(it)==t }
    }
    fun Node.otherEnd() = nodes[-kmer]!!

    fun constructGraphFromKmerSet(kmerfile: File) {
        kmerfile.forEachLine {
            val (kmer,_) = it.split('\t')
            val ckmer = Util.canonical(kmer).first
            nodes[-ckmer] = Node(-ckmer)
            nodes[ ckmer] = Node( ckmer)
        }
        if (nodes.size==0) {
            println("[Finish] No exclusive sequence. All k-mers in positive sample(+ files) occur in given negative sample(- files)")
            System.exit(1)
        }
        println("[Kmer] ${nodes.size/2} Kmers in reads containing exclusive Kmers")
        val placeholder = HashSet<Long>()
        for ((kmer,_) in nodes) if (kmer>0) {
            val prefix = kmer ushr 2
            val suffix = (kmer and (mask ushr 2)) shl 2
            for (c in 0..3){
                val (cp, _) = Util.canonical((c.toLong() shl (2*K-2))+prefix)
                val (sc, _) = Util.canonical(suffix+c)
                if (cp in nodes) nodes[kmer]!!.edges.add(Edge(placeholder, Util.decodeTable[c xor 2]))
                if (sc in nodes) nodes[-kmer]!!.edges.add(Edge(placeholder, Util.decodeTable[c]))
            }
        }
        println("</constructGraphFromKmerSet> de bruijn graph construction finished. Edge: ${nodes.map { it.value.edges.size } .sum()}")
    }
    private fun pathContraction() {
        fun compactChain(sKmer:Long): Edge{
            var kmer = sKmer
            val edge = Edge(HashSet<Long>(),StringBuilder())
            while (true){
                edge.seq.append(Util.decodeTable[(if (kmer<0) (-kmer % 4) else kmer.ushr(2*K-2).xor(2)).toInt()])
                val node = nodes[kmer]!!
                if (node.edges.size!=1 || node.otherEnd().edges.size!=1) break
                kmer = node.dest(0)
            }
            return edge
        }
        for ((_,node) in nodes)
            if (node.edges.size!=1 || node.otherEnd().edges.size!=1)
                for ( i in node.edges.indices) {
                    val t = node.dest(i)
                    node.edges[i] = compactChain(t)
                }
        for (kmer in nodes.keys.filter { it > 0 }) // remove nodes in middle of lines
            if (nodes[kmer]!!.edges.size==1 && nodes[-kmer]!!.edges.size==1) {
                nodes.remove(kmer)
                nodes.remove(-kmer)
            }
    }
    fun loadReadHopsOnGraph(file:File) {
        var readID = 0L
        val seqReader = ReadFileReader(file,"Load reads hops on De Bruijn Graph")
        while (seqReader.nextRead()){
            readID += 1
            val seq = seqReader.seq!!
            val hops = (K .. seq.length).mapNotNull {
                val node = nodes[Util.canonical(seq.substring(it-K,it)).let { if (it.second==0) it.first else -it.first}]
                if (node!=null) Pair(it,node) else null
            }.toMutableList()
            if (hops.isEmpty()) continue
            if (hops.first().first!=K){
                val e = hops.first().second.edges.find { it.seq[0]== complement[seq[hops.first().first-K-1]] }
                if (e!=null) hops.add(0,Pair(-1,nodes[hops.first().second.dest(e)]!!))
            }
            if (hops.last().first!=seq.length){
                val e = hops.last().second.otherEnd().edges.find { it.seq[0]== seq[hops.last().first] }
                if (e!=null) hops.add(Pair(-1,nodes[-hops.last().second.otherEnd().dest(e)]!!))
            }
            if (hops.size<=2) continue
            for (i in 1 until hops.size){
                val code = if (i==1) 3 else if (i==hops.size-1) 2 else 0
                hops[i].second.findEdgeTo(hops[i-1].second.kmer)?.readSet?.add((readID shl 2)+code)
                hops[i-1].second.otherEnd().findEdgeTo(hops[i].second.otherEnd().kmer)?.readSet?.add((readID shl 2)+(code xor 1))
            }
        }
    }

    fun getContigByTraverse() {
        val visit = HashMap<Long,Int>()
        fun dfs(sKmer:Long):StringBuilder {
            var kmer = sKmer
            val sb = StringBuilder()
            val readSet = HashSet<Long>()
            var readCnt = 0
            while (true){
                if (visit[Math.abs(kmer)]?:-1==readCnt) break else visit[Math.abs(kmer)] = readCnt
                val node = nodes[kmer]!!
                val edge = node.edges
                if (edge.isEmpty()) break
                val bestBranch = edge.maxBy{e-> (e.readSet.map { it ushr 2 } intersect readSet).size}!!
                readSet.addAll(bestBranch.readSet.filter { it.and(3)==2L }.map { it ushr 2 })
                val readEnd = bestBranch.readSet.filter { it.and(3)==3L }.map { it ushr 2 }
                if (bestBranch.readSet.removeIf {readSet.contains(it ushr 2)}) readCnt += 1
                readSet.removeAll(readEnd)
                sb.append(if (kmer<0) bestBranch.seq else UtilString.reverse(bestBranch.seq))
                val tKmer = node.dest(bestBranch)
                kmer = tKmer
            }
            return sb
        }
        val leafs = arrayListOf<Pair<Long,Int>>()
        for ((skmer,node) in nodes)
            if (node.edges.isEmpty()) {
                leafs.add(Pair(skmer,node.otherEnd().edges.map { it.seq.length }.max()?:0 ))
            }
        leafs.sortBy { -it.second }
        for ((skmer,_) in leafs)
            if (!visit.contains(Math.abs(skmer))){
                val s = Util.decode(if (skmer<0) Util.reverse(Math.abs(skmer)) else skmer) + dfs(-skmer).toString()
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