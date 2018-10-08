package edu.nus

import java.util.zip.GZIPInputStream


class Assembler(val settings: Desirable) {
    class Node(val kmer:String, val end:Int, var readSet: HashSet<Int>) {
        val edges:ArrayList<Edge> = arrayListOf<Edge>()
        fun dest(i:Int):Pair<String,Int>{
            val t = if (end==0) (edges[i].seq+kmer).take(kmer.length) else (kmer+edges[i].seq).takeLast(kmer.length)
            return UtilString.canonical(t)
        }
        fun dest(e:Edge):Pair<String,Int> = dest(edges.indexOf(e))
    }
    data class Edge(val readSet: HashSet<Int>, var seq:String)

    val nodes = HashMap<Pair<String,Int>,Node>()
    val contigs = mutableListOf<String>()
    fun constructGraphFromKmerSet(kmers: HashMap<Long, Long>) {
        for ((kmerCode,_) in kmers){
            val kmer = UtilString.canonical(Util.decode(kmerCode)).first
            val prefix = kmer.substring(0,kmer.length-1)
            val suffix = kmer.substring(1)
            val readset = HashSet<Int>()
            nodes[Pair(kmer,0)] = Node(kmer, 0, readset)
            nodes[Pair(kmer,1)] = Node(kmer, 1, readset)
            for (c in "ATCG"){
                val (cp, cpr) = UtilString.canonical(c+prefix)
                val (sc, scr) = UtilString.canonical(suffix+c)
                if (Util.canonical(cp) in kmers) nodes[Pair(kmer,0)]!!.edges.add(Edge(HashSet<Int>(),c.toString()))
                if (Util.canonical(sc) in kmers) nodes[Pair(kmer,1)]!!.edges.add(Edge(HashSet<Int>(),c.toString()))
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
            if (node.edges.size!=1 || nodes[Pair(kmer,end xor 1)]!!.edges.size!=1) return Edge(HashSet<Int>(),tip)
            val (tkmer,reversed) = node.dest(0)
            val edge = compactChain(tkmer,end xor reversed)
            edge.seq += tip
            node.readSet.forEach{ edge.readSet.remove(it) || edge.readSet.add(it) }
            return edge
        }
        fun dfs(kmer:String, end:Int):StringBuilder {
            visit.add(kmer)
            val node = nodes[Pair(kmer,end)]!!
            if (node.edges.isEmpty()) return StringBuilder()
            val bestBranch = node.edges.maxBy{ (readSet intersect it.readSet).size}!!
            val (tkmer,reversed) = node.dest(bestBranch)
            bestBranch.readSet.forEach{ readSet.remove(it) || readSet.add(it) }
            node.readSet.forEach{ readSet.remove(it) || readSet.add(it) }
            node.edges.remove(bestBranch)
            return dfs(tkmer,end xor reversed).append(if (end==0) bestBranch.seq else UtilString.reverse(bestBranch.seq))
        }
        for ((pair,node) in nodes)
            if (node.edges.size!=1 || nodes[Pair(pair.first,1-pair.second)]!!.edges.size!=1)
                for ( i in node.edges.indices) {
                    val (t,end) = node.dest(i)
                    node.edges[i] = compactChain(t,pair.second xor end)
                    if (pair.second==1) node.edges[i].seq = UtilString.reverse(node.edges[i].seq)
                }
        for ((pair,node) in nodes)
            if (node.edges.isEmpty() && !visit.contains(pair.first)){
                readSet.clear()
                val s = dfs(pair.first,pair.second xor 1)
                        .append(if (pair.second==0) UtilString.reverse(pair.first) else pair.first)
                contigs.add(s.toString())
            }
        println("[Assembly] ${contigs.size} contigs generated")
    }

    private fun loadReadsOnKmer() {
        val k = settings.K
        for (file in settings.files[0]){
            val seqReader = ReadFileReader(file,"Load Reads On De Bruijn Graph")
            while (seqReader.nextRead()){
                val code = (file.absolutePath+seqReader.readname).hashCode()
                val seq = seqReader.seq!!
                for (i in 0..(seq.length-k)){
                    val ckmer = UtilString.canonical(seq.substring(i,i+k)).first
                    if (Pair(ckmer,0) in nodes) {
                        nodes[Pair(ckmer,0)]!!.readSet.add(code)
                        break
                    }
                }
                for (i in (seq.length-k) downTo 0){
                    val ckmer = UtilString.canonical(seq.substring(i,i+k)).first
                    if (Pair(ckmer,0) in nodes) {
                        if (code in nodes[Pair(ckmer,0)]!!.readSet)
                            nodes[Pair(ckmer,0)]!!.readSet.remove(code)
                        else
                            nodes[Pair(ckmer,0)]!!.readSet.add(code)
                        break
                    }
                }
            }
        }
    }

    fun writeContigsToFasta(lengthLimit:Int = 0){
        println("Congtig > $lengthLimit: ${contigs.count { it.length>lengthLimit }}")
        val writer = settings.output.writer()
        for ((counter, contig) in contigs.filter { it.length>lengthLimit }.withIndex()){
            writer.write(">ExSeq${counter}_${contig.length}\n$contig\n")
        }
        writer.close()
    }
}