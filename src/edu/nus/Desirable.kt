package edu.nus

import com.github.ajalt.clikt.core.*
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.arguments.multiple
import com.github.ajalt.clikt.parameters.arguments.transformAll
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.*
import edu.gatech.kanalyze.module.KAnalyzeModule
import java.io.File

object Desirable : CliktCommand(help =
"""FILES is a list of sequencing files (fasta or fastq, may be with .gz suffix) given in following format:
    + t_1.fq t_2.fq - bg1_1.fq.gz bg1_2.fq.gz - bg_s2_1.fa bg_s2_2.fa
    where target sample is given with leading '+' and each background sample is given with leading '-'
""") {
    val threads by option(help = "Threads to use, default by the CPU core number").int().default(Runtime.getRuntime().availableProcessors())
    val K by option("-k", help = "Kmer size").int().default(31).validate { if (it>31 || it%2==0) fail("only support odd K < 32") else Util.K = it }
    val minLength by option("-m", help = "minimum contig length to report").int().default(200)
    val output by option("-o", help = "Output file path").file().default(File("DesirableOut.fasta"))
    val tmpDir: File by option("-d", help = "Temporary folder for intermediate results").file(exists = false).default(File("DesirableTmp"))
    val files: Array<ArrayList<File>> by argument().multiple().transformAll { tokenizer(it) }
    val platform = System.getProperty("os.name").split(" ")[0]
    val jarPath = File(this.javaClass.getResource("").path.split(':').last().split('!').first()).parent
    val KMCDir = when {"indows" in platform -> "win"; "ac" in platform -> "mac"; else -> "linux" }
    val KMCPath = File(jarPath,KMCDir+File.separatorChar+if (platform.contains("indows")) "kmer_counter.exe" else "kmc")
    val KMCToolsPath = File(jarPath,KMCDir+File.separatorChar+if (platform.contains("indows")) "kmc_tools.exe" else "kmc_tools")

    fun tokenizer(control: List<String>):Array<ArrayList<File>> {
        val illegalFile = control.filter { it!="+" && it!="-" }.filter { !File(it).isFile }
        if (illegalFile.isNotEmpty()) throw PrintMessage("Wrong file path:\n${illegalFile.joinToString(separator="\n")}")
        val fileGroup = Array<ArrayList<File>>(control.count { it=="+" || it=="-" },{ arrayListOf<File>() })
        if (control[0]!="+") throw PrintMessage("Sequence files for target sample should be specified before others, with leading '+' ")
        var row = -1
        for (token in control)
            if (token=="+" || token=="-") row += 1 else fileGroup[row].add(File(token))
        return fileGroup
    }
    fun logRuntime(module:String, info:String="", cmd:String){
        println("============== $module ==============\n$info")
        val runtime = kotlin.system.measureTimeMillis{
            val builder = ProcessBuilder(cmd.split(" "))
            builder.redirectOutput(ProcessBuilder.Redirect.INHERIT)
            builder.redirectError(ProcessBuilder.Redirect.INHERIT)
            val p = builder.start()
            p.waitFor()
        }
        println("[OK] finished in ${runtime/1000} seconds")
    }
    override fun run() {
        println("=== Info ===\n  Platform  : $platform\n  Thread    : $threads")
        for ((i,fileList) in files.withIndex())
            println((if (i==0) "  Target    : " else "  Background: ")+fileList.joinToString())
        tmpDir.mkdirs()
        val ko = KmerOperator(files)
        val exclusiveKmers = ko.run()
//        val exclusiveKmers = arrayListOf(File(Desirable.tmpDir,"ExtendedKmer"),File(Desirable.tmpDir,"ExclusiveRead.fastq"))
        val assembler = Assembler()
        assembler.run(exclusiveKmers[0],exclusiveKmers[1],output, minLength)
//        tmpDir.deleteRecursively()
    }
}

fun main(args: Array<String>) = Desirable.main(if (args.isEmpty()) arrayOf("--help") else args)