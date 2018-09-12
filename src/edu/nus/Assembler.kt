package edu.nus

import kotlin.coroutines.experimental.EmptyCoroutineContext.plus

object Util{
    val supl = mapOf<Char, Char>('A' to 'T', 'C' to 'G', 'G' to 'C', 'T' to 'A')
    fun reverse(kmer: String) = kmer.reversed().map { supl[it] }.joinToString(separator = "")
    fun canonical(kmer: String) = if (kmer[kmer.length/2] in "AC") Pair(kmer,0) else Pair(reverse(kmer),1)
}

typealias Node = Pair<String,Int>
class Assembler(val settings: Desirable) {
    val edge = HashMap<Node,MutableList<Node>>()
    val contigs = mutableListOf<String>()
    fun constructGraphFromKmerSet(rawKmers: HashSet<String>) {
        val kmers = rawKmers.map { Util.canonical(it).first } .toSet()
        for (kmer in kmers){
            val prefix = kmer.substring(0,kmer.length-1)
            val suffix = kmer.substring(1)
            edge[Pair(kmer,0)] = mutableListOf<Node>()
            edge[Pair(kmer,1)] = mutableListOf<Node>()
            for (c in "ATCG"){
                val (cp,cpEnd) = Util.canonical(c+prefix)
                val (sc,scEnd) = Util.canonical(suffix+c)
                if (cp in kmers) edge[Pair(kmer,0)]!!.add(Pair(cp,1-cpEnd))
                if (sc in kmers) edge[Pair(kmer,1)]!!.add(Pair(sc,scEnd))
            }
        }
        println("[Assembly] de bruijn graph construction finished\nEdge: ${edge.map { it.value.size } .sum()}")
        for (kmer in kmers){
            if (edge[Pair(kmer,0)]!!.size>1 || edge[Pair(kmer,1)]!!.size>1){
                edge[Pair(kmer,0)] = mutableListOf<Node>()
                edge[Pair(kmer,1)] = mutableListOf<Node>()
            }
        }
        println("[Assembly] Simple chain check finished\nEdge: ${edge.map { it.value.size } .sum()}")
    }
    fun getContigByTraverse() {
        val visit = HashSet<String>()
        fun dfs(node:Node):StringBuilder {
            val (kmer, end) = node
            visit.add(node.first)
            if (edge[Pair(kmer,1-end)]!!.isEmpty())
                return StringBuilder(if (end==0) Util.reverse(kmer) else kmer)
            return dfs(edge[Pair(kmer,1-end)]!![0]).append(if (end==0) Util.supl[kmer[0]] else kmer.last())
        }
        for ((node,neighbor) in edge)
            if (neighbor.isEmpty() && !visit.contains(node.first))
                contigs.add(dfs(node).toString())
        println("[Assembly] ${contigs.size} contigs generated")
    }
    fun writeContigsToFasta(lengthLimit:Int = 0){
        val writer = settings.output.writer()
        for ((counter, contig) in contigs.filter { it.length>lengthLimit }.withIndex()){
            writer.write(">ExSeq${counter}_${contig.length}\n$contig\n")
        }
        writer.close()
    }
}