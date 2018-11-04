package edu.nus

import com.github.ajalt.clikt.core.*
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.arguments.multiple
import com.github.ajalt.clikt.parameters.arguments.transformAll
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.*
import edu.gatech.kanalyze.module.KAnalyzeModule
import java.io.File

class Desirable : CliktCommand(help =
"""FILES is a list of sequencing files (fasta or fastq, may be with .gz suffix) given in following format:
    + t_1.fq t_2.fq - bg1_1.fq.gz bg1_2.fq.gz - bg_s2_1.fa bg_s2_2.fa
    where target sample is given with leading '+' and each background sample is given with leading '-'
""") {
    val threads by option(help = "Threads to use, default by the CPU core number").int().default(Runtime.getRuntime().availableProcessors())
    val K by option("-k", help = "Kmer size").int().default(31).validate { if (it>32) fail("only support K <= 32") else Util.K = it }
    val minLength by option("-m", help = "minimum contig length to report").int().default(200)
    val output by option("-o", help = "Output file path").file().default(File("DesirableOut.fasta"))
    val tmpDir: File by option("-d", help = "Temporary folder for intermediate results").file(exists = false).default(File("DesirableTmp"))
    val files: Array<ArrayList<File>> by argument().multiple().transformAll { tokenizer(it) }
    var countFiles: Array<File> = arrayOf()

    fun tokenizer(control: List<String>):Array<ArrayList<File>> {
        val illegalFile = control.filter { it!="+" && it!="-" }.filter { !File(it).isFile }
        if (illegalFile.isNotEmpty()) throw PrintMessage("Wrong file path:\n${illegalFile.joinToString(separator="\n")}")
        val fileGroup = Array<ArrayList<File>>(control.count { it=="+" || it=="-" },{ arrayListOf<File>() })
        if (control[0]!="+") throw PrintMessage("Sequence files for target sample should be specified before others, with leading '+' ")
        var row = -1
        for (token in control)
            if (token=="+" || token=="-") row += 1 else fileGroup[row].add(File(token))
        for ((i,fileList) in fileGroup.withIndex())
            println((if (i==0) "Target     : " else "Background : ")+fileList.joinToString())
        return fileGroup
    }
    fun countKmerInFile(files: List<File>):File {
        val filename = files.first().name
        val countFile = File(tmpDir, "$filename.kc")
        val format = (if (filename.contains("fastq") || filename.contains("fq")) "fastq" else "fasta") + if (filename.endsWith(".gz")) "gz" else ""
        val arguments = "count -rcanonical " +
                "-t $threads " +
                "-k $K " +
                "-c kmercount:2 " +
                "-o ${countFile.path} " +
                "-f $format ${files.joinToString(separator=" ") { it.path }}"
        println("[KmerCount] ====== start ======\n[KmerCount] Command: $arguments")
        val runtime = kotlin.system.measureTimeMillis { KAnalyzeModule.main(arguments.split(" ").toTypedArray()) }
        println("[KmerCount] finished in ${runtime/1000} seconds")
        return countFile
    }

    override fun run() {
        tmpDir.mkdirs()
        countFiles = files.map { countKmerInFile(it) }.toTypedArray()
        val ka = KmerAnalyzer(this)
        val exclusiveKmers = ka.filterTargetByControl()
        val assembler = Assembler(this)
        assembler.constructGraphFromKmerSet(exclusiveKmers)
        assembler.getContigByTraverse()
        assembler.writeContigsToFasta(minLength)
//        tmpDir.deleteRecursively()
    }
}

fun main(args: Array<String>) = Desirable().main(if (args.isEmpty()) arrayOf("--help") else args)