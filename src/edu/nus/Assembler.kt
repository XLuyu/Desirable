package edu.nus

import java.io.File


class Assembler() {
    val nodes = HashMap<Pair<String,Int>,Node>()
    val contigs = mutableListOf<String>()

    data class Edge(val readStart: HashSet<Int>, val readSet: HashSet<Int>, val readEnd: HashSet<Int>, var seq:StringBuilder){
        // seq is store as if append after this kmer
        constructor(s: HashSet<Int>, t: HashSet<Int>, e: HashSet<Int>, c:Char) : this(s,t,e,StringBuilder(c.toString()))
    }
    class Node(val kmer:String, val end:Int, val readStart: HashSet<Int>, val readEnd: HashSet<Int>) {
        val edges:ArrayList<Edge> = arrayListOf<Edge>()
        fun dest(e:Edge):Pair<String,Int> {
            val t = ((if (end==0) kmer.revCompl() else kmer)+e.seq).takeLast(kmer.length)
            return UtilString.canonical(t).let { Pair(it.first,it.second xor 1) }
        }
        fun dest(i:Int):Pair<String,Int>{
            val t = ((if (end==0) kmer.revCompl() else kmer)+edges[i].seq).takeLast(kmer.length)
            return UtilString.canonical(t).let { Pair(it.first,it.second xor 1) }
        }
    }
    fun Node.otherEnd() = nodes[Pair(kmer, end xor 1)]!!

    fun constructGraphFromKmerSet(kmerfile: File) {
        val kmers = HashMap<Long, Long>()
        kmerfile.forEachLine {
            val (kmer,count) = it.split('\t')
            kmers[Util.canonical(kmer)] = count.toLong()
        }
        println("[Kmer] ${kmers.size} Kmers in reads containing exclusive Kmers")
        for ((kmerCode,_) in kmers) {
            val kmer = UtilString.canonical(Util.decode(kmerCode)).first
            nodes[Pair(kmer, 0)] = Node(kmer, 0, HashSet<Int>(), HashSet<Int>())
            nodes[Pair(kmer, 1)] = Node(kmer, 1, HashSet<Int>(), HashSet<Int>())
        }
        for ((kmerCode,_) in kmers) {
            val kmer = UtilString.canonical(Util.decode(kmerCode)).first
            val prefix = kmer.substring(0, kmer.length - 1)
            val suffix = kmer.substring(1)
            for (c in "ATCG"){
                val (cp, cpr) = UtilString.canonical(c+prefix)
                val (sc, scr) = UtilString.canonical(suffix+c)
                if (Util.canonical(cp) in kmers && Util.canonical(cp)<kmerCode) {
                    val newStart = HashSet<Int>()
                    val newEnd = HashSet<Int>()
                    val newSet = HashSet<Int>()
                    nodes[Pair(kmer, 0)]!!.edges.add(Edge(newStart,newSet,newEnd, complement[c]!!))
                    nodes[Pair(cp, 1-cpr)]!!.edges.add(Edge(newStart,newSet,newEnd, kmer.last()))
                }
                if (Util.canonical(sc) in kmers && Util.canonical(sc)<kmerCode) {
                    val newStart = HashSet<Int>()
                    val newEnd = HashSet<Int>()
                    val newSet = HashSet<Int>()
                    nodes[Pair(kmer,1)]!!.edges.add(Edge(newStart,newSet,newEnd, c))
                    nodes[Pair(sc, scr)]!!.edges.add(Edge(newStart,newSet,newEnd, complement[kmer[0]]!!))
                }
            }
        }
        println("[Assembly] de bruijn graph construction finished\nEdge: ${nodes.map { it.value.edges.size } .sum()}")
    }

    fun loadReadEndsOnGraph(file:File) {
        val k = Desirable.K
        val seqReader = ReadFileReader(file,"Load reads ends on De Bruijn Graph")
        while (seqReader.nextRead()){
            val code = (file.absolutePath+seqReader.readname).hashCode()
            val seq = seqReader.seq!!
            for (i in 0..(seq.length-k)){
                val node = nodes[UtilString.canonical(seq.substring(i,i+k))]
                if (node!=null) {
                    node.readEnd.add(code)
                    node.otherEnd().readStart.add(code)
                    break
                }
            }
            for (i in (seq.length-k) downTo 0){
                val node = nodes[UtilString.canonical(seq.substring(i,i+k))]
                if (node!=null) {
                    node.readStart.add(code)
                    node.otherEnd().readEnd.add(code)
                    break
                }
            }
        }
    }

    private fun pathContraction() {
        fun compactChain(sKmer:String, sEnd:Int): Edge{
            var kmer = sKmer
            var end = sEnd
            val edge = Edge(HashSet<Int>(),HashSet<Int>(),HashSet<Int>(),StringBuilder())
            while (true){
                val node = nodes[Pair(kmer,end)]!!
                edge.seq.append(if (end==1) kmer.last() else complement[kmer[0]])
                if (node.edges.size!=1 || node.otherEnd().edges.size!=1) break
                edge.readStart.addAll(node.readStart)
                edge.readEnd.addAll(node.readEnd)
                val (tKmer,tEnd) = node.dest(0)
                kmer = tKmer
                end = tEnd
            }
            return edge
        }
        for ((pair,node) in nodes)
            if (node.edges.size!=1 || node.otherEnd().edges.size!=1)
                for ( i in node.edges.indices) {
                    val (t,end) = node.dest(i)
                    node.edges[i] = compactChain(t,end)
                    val inside = node.edges[i].readStart.filter { node.edges[i].readEnd.remove(it) }
                    node.edges[i].readStart.removeAll(inside)
                }
        for ((kmer,end) in nodes.keys.filter { it.second==0 })
            if (nodes[Pair(kmer,0)]!!.edges.size==1 && nodes[Pair(kmer,1)]!!.edges.size==1) {
                nodes.remove(Pair(kmer,0))
                nodes.remove(Pair(kmer,1))
            }
    }
    fun loadReadHopsOnGraph(file:File) {
        val k = Desirable.K
        val seqReader = ReadFileReader(file,"Load reads hops on De Bruijn Graph")
        while (seqReader.nextRead()){
            val code = (file.absolutePath+seqReader.readname).hashCode()
            val seq = seqReader.seq!!
            for (i in k .. seq.length){
                val node = nodes[UtilString.canonical(seq.substring(i-k,i))]
                if (node!=null) {
                    if (k<i){
                        val c = complement[seq[i-k-1]]
                        node.edges.find { it.seq[0]== c }?.readSet?.add(code)
                    }
                    if (i<seq.length){
                        val c = seq[i]
                        node.otherEnd().edges.find { it.seq[0]==c }?.readSet?.add(code)
                    }
                }
            }
        }
    }

    fun getContigByTraverse() {
        val visit = HashMap<String,Int>()
        fun dfs(sKmer:String, sEnd:Int):StringBuilder {
            var kmer = sKmer
            var end = sEnd
            val sb = StringBuilder()
            val readSet = HashSet<Int>()
            var readCnt = 0
            while (true){
                if (visit[kmer]?:-1==readCnt) break else visit[kmer] = readCnt
                visit[kmer] = readCnt
                val node = nodes[Pair(kmer,end)]!!
                val edge = node.edges
                if (edge.isEmpty()) break
                readSet.addAll(node.readStart)
                readSet.removeAll(node.readEnd)
                val bestBranch = edge.maxBy{ (readSet intersect it.readSet).size}!!
                if (bestBranch.readSet.removeAll(readSet)) readCnt += 1
                readSet.addAll(bestBranch.readStart)
                readSet.removeAll(bestBranch.readEnd)
                sb.append(if (end==1) bestBranch.seq else UtilString.reverse(bestBranch.seq))
                val (tKmer,tEnd) = node.dest(bestBranch)
                kmer = tKmer
                end = tEnd
            }
            return sb
        }
        val leafs = arrayListOf<Pair<Pair<String,Int>,Int>>()
        for ((pair,node) in nodes)
            if (node.edges.isEmpty()) {
                leafs.add(Pair(pair,node.otherEnd().edges.map { it.seq.length }.max()?:0 ))
            }
        leafs.sortBy { -it.second }
        for ((pair,_) in leafs)
            if (!visit.contains(pair.first)){
                val s = (if (pair.second==1) pair.first.revCompl() else pair.first) + dfs(pair.first,pair.second xor 1).toString()
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
        loadReadEndsOnGraph(reads)
        pathContraction()
        loadReadHopsOnGraph(reads)
        getContigByTraverse()
        writeContigsToFasta(out,minLength)
    }
}