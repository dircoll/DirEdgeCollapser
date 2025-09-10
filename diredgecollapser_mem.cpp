#include "DIR_EDGE_COLLAPSER_MEM.H" // Update this to match your header filename
#include <fstream>
#include <iostream>
#include <sstream>
#include <tuple>
#include <vector>
#include <set>
#include <chrono>
#include <iomanip>
#include <cstring>

// Platform-specific includes for memory tracking
#ifdef _WIN32
    #include <windows.h>
    #include <psapi.h>
#elif __linux__
    #include <sys/resource.h>
    #include <unistd.h>
#elif __APPLE__
    #include <mach/mach.h>
#endif

// Memory tracking functions
class MemoryTracker {
private:
    size_t initial_memory;
    size_t peak_memory;
    
public:
    MemoryTracker() : initial_memory(0), peak_memory(0) {
        initial_memory = getCurrentMemoryUsage();
        peak_memory = initial_memory;
    }
    
    // Get current memory usage in bytes
    static size_t getCurrentMemoryUsage() {
#ifdef _WIN32
        PROCESS_MEMORY_COUNTERS_EX pmc;
        GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
        return pmc.WorkingSetSize;
#elif __linux__
        long rss = 0L;
        FILE* fp = fopen("/proc/self/status", "r");
        if (fp != NULL) {
            char line[128];
            while (fgets(line, 128, fp) != NULL) {
                if (strncmp(line, "VmRSS:", 6) == 0) {
                    char* pEnd;
                    rss = strtol(line + 6, &pEnd, 10) * 1024; // Convert KB to bytes
                    break;
                }
            }
            fclose(fp);
        }
        return rss;
#elif __APPLE__
        struct mach_task_basic_info info;
        mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
        if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, 
                     (task_info_t)&info, &infoCount) == KERN_SUCCESS) {
            return info.resident_size;
        }
        return 0;
#else
        return 0; // Unsupported platform
#endif
    }
    
    void updatePeak() {
        size_t current = getCurrentMemoryUsage();
        if (current > peak_memory) {
            peak_memory = current;
        }
    }
    
    size_t getInitialMemory() const { return initial_memory; }
    size_t getPeakMemory() const { return peak_memory; }
    size_t getCurrentIncrease() const { 
        return getCurrentMemoryUsage() - initial_memory; 
    }
    size_t getPeakIncrease() const { 
        return peak_memory - initial_memory; 
    }
};

// Helper function to format memory size
std::string formatMemory(size_t bytes) {
    const char* units[] = {"B", "KB", "MB", "GB"};
    int unit_index = 0;
    double size = static_cast<double>(bytes);
    
    while (size >= 1024.0 && unit_index < 3) {
        size /= 1024.0;
        unit_index++;
    }
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << size << " " << units[unit_index];
    return oss.str();
}

// Define a pair of source and destination for easy comparison
using Vrtx_Pair = std::pair<Vertex, Vertex>;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
        return 1;
    }

    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];

    // Initialize memory tracking
    MemoryTracker memTracker;
    
    std::cout << "Initial memory usage: " << formatMemory(memTracker.getInitialMemory()) << std::endl;

    std::ifstream inputFile(inputFileName);
    std::ofstream outputFile(outputFileName);

    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open input file " << inputFileName << std::endl;
        return 1;
    }

    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputFileName << std::endl;
        return 1;
    }

    std::string line;
    std::vector<std::string> dim0_lines;
    std::vector<Edge> edges, c_int_edges;
    std::set<Vrtx_Pair> edgeSet;
    bool in_dim0 = false, in_dim1 = false;

    // Parse the input file
    std::cout << "Reading input file..." << std::endl;
    while (std::getline(inputFile, line)) {
        if (line.find("dim 0") != std::string::npos) {
            dim0_lines.push_back(line);
            in_dim0 = true;
            in_dim1 = false;
            continue;
        } else if (line.find("dim 1") != std::string::npos) {
            dim0_lines.push_back(line); // Also save the "dim 1:" line
            in_dim0 = false;
            in_dim1 = true;
            continue;
        }

        if (in_dim0 && !line.empty()) {
            dim0_lines.push_back(line);
        }

        if (in_dim1 && !line.empty()) {
            std::stringstream ss(line);
            int src, dest;
            double weight;
            if (ss >> src >> dest >> weight) {
                Vrtx_Pair uv = std::make_pair(src, dest);
                Vrtx_Pair vu = std::make_pair(dest, src);
                
                // Only add edge if reverse edge not already present
                if (edgeSet.find(vu) == edgeSet.end()) {
                    edges.emplace_back(src, dest, weight);
                    edgeSet.insert(uv);
                }
            }
        }
    }

    inputFile.close();
    
    // Update memory after loading
    memTracker.updatePeak();
    std::cout << "Memory after loading graph: " << formatMemory(MemoryTracker::getCurrentMemoryUsage()) << std::endl;
    std::cout << "Memory used for loading: " << formatMemory(memTracker.getCurrentIncrease()) << std::endl;

    std::cout << "\nTotal number of initial edges: " << edges.size() << std::endl;
    
    // Time the algorithm execution
    std::cout << "\nRunning directed flag complex edge collapser algorithm..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Track memory before algorithm
    size_t memory_before_algo = MemoryTracker::getCurrentMemoryUsage();
    
    // Run the directed flag complex edge collapser algorithm
    {
        CoreDirFlagFiltration cff(edges);
        
        // Update peak memory after graph construction
        memTracker.updatePeak();
        std::cout << "Memory after graph construction: " << formatMemory(MemoryTracker::getCurrentMemoryUsage()) << std::endl;
        
        cff.process_edges(c_int_edges);
        
        // Update peak memory after processing
        memTracker.updatePeak();
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // Memory statistics
    size_t memory_after_algo = MemoryTracker::getCurrentMemoryUsage();
    size_t algo_memory_usage = memory_after_algo - memory_before_algo;
    
    std::cout << "\nAlgorithm execution completed!" << std::endl;
    std::cout << "================== Performance Statistics ==================" << std::endl;
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
    std::cout << "Initial edges: " << edges.size() << std::endl;
    std::cout << "Filtered edges: " << c_int_edges.size() << std::endl;
    
    // Calculate reduction percentage
    double reduction = (1.0 - static_cast<double>(c_int_edges.size()) / edges.size()) * 100.0;
    std::cout << "Edge reduction: " << std::fixed << std::setprecision(2) << reduction << "%" << std::endl;
    
    std::cout << "\n================== Memory Usage Statistics =================" << std::endl;
    std::cout << "Initial memory: " << formatMemory(memTracker.getInitialMemory()) << std::endl;
    std::cout << "Peak memory: " << formatMemory(memTracker.getPeakMemory()) << std::endl;
    std::cout << "Peak memory increase: " << formatMemory(memTracker.getPeakIncrease()) << std::endl;
    std::cout << "Current memory: " << formatMemory(memory_after_algo) << std::endl;
    std::cout << "Algorithm memory usage: " << formatMemory(algo_memory_usage) << std::endl;
    
    // Memory efficiency metrics
    double bytes_per_edge = static_cast<double>(memTracker.getPeakIncrease()) / edges.size();
    std::cout << "Memory per edge: " << std::fixed << std::setprecision(2) << bytes_per_edge << " bytes/edge" << std::endl;
    std::cout << "==========================================================" << std::endl;

    // Write output file
    std::cout << "\nWriting output file..." << std::endl;
    
    // Write dim 0 section first (vertices)
    bool dim1_written = false;
    for (const auto& l : dim0_lines) {
        if (l.find("dim 1") != std::string::npos) {
            dim1_written = true;
        }
        outputFile << l << std::endl;
    }

    // Write dim 1 section header if not already written
    if (!dim1_written) {
        outputFile << "dim 1:" << std::endl;
    }

    // Write the filtered edges
    for (const auto& edge : c_int_edges) {
        outputFile << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << std::endl;
    }

    outputFile.close();

    std::cout << "\nFile successfully processed and written to " << outputFileName << std::endl;
    
    // Final memory summary
    std::cout << "\nFinal memory usage: " << formatMemory(MemoryTracker::getCurrentMemoryUsage()) << std::endl;
    std::cout << "Total memory released: " << formatMemory(memory_before_algo - MemoryTracker::getCurrentMemoryUsage()) << std::endl;
    
    return 0;
}
