package edu.nus

import java.io.File
import java.util.zip.GZIPInputStream


class Assembler() {
    val nodes = HashMap<Pair<String,Int>,Node>()
    val contigs = mutableListOf<String>()

    class Node(val kmer:String, val end:Int, val readStart: HashSet<Int>, val readEnd: HashSet<Int>) {
        val edges:ArrayList<Edge> = arrayListOf<Edge>()
        fun dest(i:Int):Pair<String,Int>{
            val t = if (end==0) (edges[i].seq+kmer).take(kmer.length) else (kmer+edges[i].seq).takeLast(kmer.length)
            return UtilString.canonical(t)
        }
        fun dest(e:Edge):Pair<String,Int> = dest(edges.indexOf(e))
    }
    data class Edge(val readSet: HashSet<Int>, var seq:String)

    fun constructGraphFromKmerSet(kmers: HashMap<Long, Long>) {
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
                    val newSet = HashSet<Int>()
                    nodes[Pair(kmer, 0)]!!.edges.add(Edge(newSet, c.toString()))
                    nodes[Pair(cp, 1-cpr)]!!.edges.add(Edge(newSet, if (cpr==0) kmer.last().toString() else UtilString.supl[kmer.last()].toString()))
                }
                if (Util.canonical(sc) in kmers && Util.canonical(sc)<kmerCode) {
                    val newSet = HashSet<Int>()
                    nodes[Pair(kmer,1)]!!.edges.add(Edge(newSet,c.toString()))
                    nodes[Pair(sc, scr)]!!.edges.add(Edge(newSet, if (scr==0) kmer.first().toString() else UtilString.supl[kmer.first()].toString()))
                }
            }
        }
        println("[Assembly] de bruijn graph construction finished\nEdge: ${nodes.map { it.value.edges.size } .sum()}")
    }
    var debug = false
    fun getContigByTraverse() {
        val visit = HashSet<String>()
        val readSet = HashSet<Int>()
        fun compactChain(kmer:String, end:Int):Edge{
            val tip = (if (end==0) kmer[0] else UtilString.supl[kmer.last()]).toString()
            val node = nodes[Pair(kmer,end)]!!
            if (node.edges.size!=1) return Edge(HashSet<Int>(),tip)
            val (tkmer,reversed) = node.dest(0)
            val edge = compactChain(tkmer,end xor reversed)
            edge.seq += tip
            return edge
        }
        fun dfs(kmer:String, end:Int):StringBuilder {
            visit.add(kmer)
            val node = nodes[Pair(kmer,end)]!!
            val edge = node.edges
            if (edge.isEmpty()) return StringBuilder()
            readSet.addAll(node.readStart)
            readSet.removeAll(node.readEnd)
            val bestBranch = edge.maxBy{ (readSet intersect it.readSet).size}!!
            val (tkmer,reversed) = node.dest(bestBranch)
            bestBranch.readSet.removeAll(readSet)
            return dfs(tkmer,end xor reversed).append(if (end==0) bestBranch.seq else UtilString.reverse(bestBranch.seq))
        }
        val leafs = arrayListOf<Pair<Pair<String,Int>,Int>>()
//        for ((pair,node) in nodes)
//            if (node.edges.size!=1 || nodes[Pair(pair.first,1-pair.second)]!!.edges.size!=1) {
//                for (i in node.edges.indices) {
//                    val (t, end) = node.dest(i)
//                    node.edges[i] = compactChain(t, pair.second xor end)
//                    if (pair.second == 1) node.edges[i].seq = UtilString.reverse(node.edges[i].seq)
//                }
//            }
        for ((pair,node) in nodes)
            if (node.edges.isEmpty()) {
                leafs.add(Pair(pair,compactChain(pair.first,pair.second xor 1).seq.length-1))
            }
        leafs.sortBy { -it.second }
        for ((pair,dist) in leafs)
            if (!visit.contains(pair.first)){
                readSet.clear()
                val s = dfs(pair.first,pair.second xor 1)
                        .append(if (pair.second==0) UtilString.reverse(pair.first) else pair.first)
                contigs.add(s.toString())
            }
        println("[Assembly] ${contigs.size} contigs generated")
    }

    fun loadReadsOnKmer(files:Array<File>) {
        val k = Desirable.K
        for (file in files){
            val seqReader = ReadFileReader(file,"Load Reads On De Bruijn Graph")
            while (seqReader.nextRead()){
                val code = (file.absolutePath+seqReader.readname).hashCode()
                val seq = seqReader.seq!!
                for (i in 0..(seq.length-k)){
                    val node = UtilString.canonical(seq.substring(i,i+k))
                    if (node in nodes) {
                        nodes[node]!!.readEnd.add(code)
                        nodes[Pair(node.first,1-node.second)]!!.readStart.add(code)
                        break
                    }
                }
                for (i in (seq.length-k) downTo 0){
                    val node  = UtilString.canonical(seq.substring(i,i+k))
                    if (node in nodes) {
                        nodes[node]!!.readStart.add(code)
                        nodes[Pair(node.first,1-node.second)]!!.readEnd.add(code)
                        break
                    }
                }
                var last = UtilString.canonical(seq.substring(0,k))
                for (i in k+1 .. seq.length){
                    val kmer = UtilString.canonical(seq.substring(i-k,i))
                    if (kmer in nodes && last in nodes) {
                        val node = nodes[kmer]!!
                        node.edges.find { val (mer,r) = node.dest(it); mer==last.first && kmer.second xor r==last.second }!!.readSet.add(code)
                    }
                    last = kmer
                }
            }
        }
    }

    fun writeContigsToFasta(file: File, lengthLimit:Int = 0){
        println("Congtig > $lengthLimit: ${contigs.count { it.length>lengthLimit }}")
        val writer = file.writer()
        for ((counter, contig) in contigs.filter { it.length>lengthLimit }.withIndex()){
            writer.write(">ExSeq${counter}_${contig.length}\n$contig\n")
        }
        writer.close()
    }

    fun constructGraphFromKmerSet(kmerfile: File) {
        val kmers = HashMap<Long, Long>()
        kmerfile.forEachLine {
            val (kmer,count) = it.split('\t')
            kmers[Util.canonical(kmer)] = count.toLong()
        }
        println("[Kmer] ${kmers.size} Kmers in reads containing exclusive Kmers")
        constructGraphFromKmerSet(kmers)
    }
}