# Unistra HPC Resources (CAIUS)

Created: February 12, 2026 9:19 AM
Tags: HPC
Author: Tanul Gupta
Purpose: Quick, reliable reference for choosing Slurm partition + resources so jobs schedule smoothly.

## **Unistra HPC (CAIUS)**

### **CPU Node Configurations (as of 10-06-2024)**

**Cluster totals:** 541 compute nodes • 13,078 CPU cores • ~68 TB RAM

**Rule of thumb:** ~2–8 GB RAM per core (varies by node type)

| **Nodes** | **CPU Generation** | **Architecture** | **Sockets × Cores** | **Total Cores / Node** | **RAM per Node (GB)** | **Typical RAM/Core** |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | AMD EPYC Rome Zen2 7002 | AMD | 2 × 24 | **48** | 128 / 256 | 2.7–5.3 GB |
| 63 | Intel Cascade Lake | Intel | 2 × 16 | **32** | 128 / 384 | 4–12 GB |
| 3 | Intel Cascade Lake | Intel | 2 × 8 | **16** | 128 | 8 GB |
| 127 | Intel Skylake | Intel | 2 × 12 | **24** | 96 / 192 | 4–8 GB |
| 73 | Intel Broadwell | Intel | 2 × 14 | **28** | 128 / 256 | 4.6–9.1 GB |
| 35 | Intel Haswell | Intel | 2 × 8 | **16** | 64 / 128 | 4–8 GB |
| 95 | Intel Ivy Bridge | Intel | 2 × 8 | **16** | 64 / 128 | 4–8 GB |
| 166 | Intel Sandy Bridge | Intel | 2 × 8 | **16** | 64 | 4 GB |
| 4 | Intel Sandy Bridge | Intel | 2 × 6 | **12** | 32 | 2.7 GB |
| 16 | Intel Xeon Nehalem | Intel | 2 × 6 | **12**  | 24 | 2 GB |
| **592** |  |  |  |  |  |  |
|  |  |  |  |  |  |  |

## **2. Main User-Facing Partitions**

### **Status Summary**

| **Partition** | **Type** | **Time Limit** | **Total Nodes** | **Reserved (Planned)** | **Unavailable** | RAM | **Typical Use** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **public** | CPU | 1 day | **347** | 28 | 4 |  | largest pool → fastest scheduling, heterogeneous hardware. |
| **grant** | CPU | 4 days | **234** | 28 | 4 | 39 TB | same machines but different accounting and longer allowed runtime. |
| **publicgpu** | GPU | 1 day | **37** | 0 | 1 |  | shared GPU queue with mixed usage. |
| **grantgpu** | GPU | 4 days | **12** | 0 | 0 |  | smaller, project-controlled GPU allocation. |

### **Practical decision flow**

Start →

- Need **>24h**? → **grant** (max 4d)
- Need **>4 nodes**? → **grant** (max 32 nodes/job)
- Doing **final long stat / production**? → prefer **grant**
- Want **fast turnaround**? → **public** (default)

### **Private / Thematic Partitions (Restricted Use)**

These partitions correspond to **specific research groups or funding lines**:

| **Partition** | **Meaning** | Time Limit |
| --- | --- | --- |
| **pri2025, pri2020, pri2016, etc.** | Year- or project-based private allocations | Infinite |
| **a2s, infochem2018** | Thematic projects | Infinite |
| **igbmcgpu, priipgsgpu, ceremagpu** | Lab-specific GPU pools | Infinite |
| **priamd2021** | Dedicated AMD nodes | Infinite |
| **pri2018gpu, pri2021gpu** | Older GPU allocations | 10-20 days |

**NOTE:** You can only use these if your account is mapped to them.

---

## **3. `public` vs `grant` (CPU)**

| **Property** | **public** | **grant** | **What It Means Practically** |
| --- | --- | --- | --- |
| Default Partition | **YES** | NO | Jobs go to public if you don’t specify -p |
| Default Time Limit | **1 day** | **4 days** | grant allows longer simulations |
| Maximum Time | 1 day | 4 days | Hard runtime ceiling |
| Max Nodes per Job | **4** | **32** | grant supports large MPI jobs |
| Total Nodes | **347** | **234** | public has more machines overall |
| Total CPUs | **9388 cores** | **6000 cores** | Larger raw capacity in public |
| Total Memory (TRES) | ~57 TB | ~39 TB | More aggregate RAM in public |
| Node Pool | Superset (includes extra nodes) | Subset of public hardware | grant uses mostly same machines |
| Oversubscription | NO | NO | Nodes are not oversubscribed; each requested core is dedicated and not shared with other jobs |
| Preemption Mode | REQUEUE | REQUEUE | Job stops → returns to queue → restarts later |
| Memory Limits | Unlimited (must request explicitly) | Unlimited (must request explicitly) | Always set --mem yourself |
| QoS | N/A | N/A | No special QoS differentiation |
| Access Restrictions | None | None | Access controlled by accounting |
| State | UP | UP | Both active partitions |

<aside>
<img src="notion://custom_emoji/2836a877-8186-4f74-8c80-16dff2ee7c78/21d16d38-5ead-80b9-bbc0-007aa06af60f" alt="notion://custom_emoji/2836a877-8186-4f74-8c80-16dff2ee7c78/21d16d38-5ead-80b9-bbc0-007aa06af60f" width="40px" />

**Structural Relationship Between Them**

- The grant partition is essentially: $\textbf{grant} \subset \textbf{public}$
- public includes **additional nodes** (e.g., 519–520, 855–861, 883–886, etc.)
- grant is a **restricted subset intended for longer production workloads**
</aside>

## Nodewise Summary

- **SLURM COMMANDS**
    
    
    **A) Nodewise list (one line per node)**
    
    ```
    sinfo -p grant -N -h -o "%N %c %m %f"
    ```
    
    - %N node name
    - %c cores
    - %m memory (MB)
    - %f features (CPU generation, etc.)
    
    **B) Nodewise summary (grouped by node type)**
    
    Count how many nodes of each hardware type exist
    
    ```
    sinfo -p grant -N -h -o "%c|%m|%f" | sort | uniq -c
    ```
    
    **C) Summary with node ranges + automatic node counts (most polished)**
    This produces a compact table with **NodeRange, #Nodes, cores, RAM(GB), features**:
    
    ```bash
    sinfo -p grant -o "%N|%c|%m|%f" | tail -n +2 | \
    while IFS='|' read -r nodelist c m f; do
      n=$(scontrol show hostnames "$nodelist" | wc -l)
      printf "%-70s %4d  %2s cores  %7.1f GB  %s\n" "$nodelist" "$n" "$c" "$(echo "$m/1024" | bc -l)" "$f"
    done
    ```
    

### 1.  `grant` partition (12 Feb 2026)

## **Partition — Node Range Summary**

| **Node Range** | **# Nodes** | **CPU Generation** | **Cores / Node** | **Memory / Node (GB)** | **Notes** |
| --- | --- | --- | --- | --- | --- |
| hpc-n[675-700,702-714,719-738,755,757-772,774-785,787-788,790,792-795,797-799,801-805,807-824,826-832,834,836-849,851,853] | 145 | Intel Skylake | 24 | 91.8 | Majority nodes (standard tier) |
| hpc-n[651-652,654-655,657-674] | 22 | Intel Skylake | 24 | 187.3 | Medium-memory Skylake |
| hpc-n593 | 1 | Intel Skylake | 24 | 91.8 | Single-node variant |
| hpc-n[523-546] | 24 | Intel Broadwell | 28 | 125.0 | Older generation |
| hpc-n[547-548,550-554,556-560] | 12 | Intel Broadwell | 28 | 250.7 | High-memory Broadwell |
| hpc-n876 | 1 | Intel Cascade Lake | 32 | 184.6 | Single newer Intel node |
| hpc-n[935-942,944-951,953-958,960-966] | 29 | AMD EPYC (epyc4) | 32 | **500.0** | Large-memory nodes |

$$
145 + 22 + 1+ 24 + 12 + 1 + 29 = 234 \text{ nodes}
$$

### 2.  `public` partition (12 Feb 2026)

## **Partition — Node Range Summary**

| **Node Range** | **#Nodes** | **CPU Generation** | **Cores / Node** | **Memory / Node (GB)** | **Notes** |
| --- | --- | --- | --- | --- | --- |
| hpc-n[519-520,523-546,561-567,569-575,578-590] | 53 | Intel Broadwell | 28 | 128 | AVX2 |
| hpc-n[547-548,550-554,556-560] | 12 | Intel Broadwell | 28 | 256 | High-memory Broadwell |
| hpc-n[591-593,631-650,675-700,702-714,719-738,755,757-772,774-785,787-788,790,792-795,797-799,801-805,807-824,826-832,834,836-849,851,853] | 168 | Intel Skylake | 24 | 96 | AVX-512 |
| hpc-n[623-628,630,651-652,654-655,657-674] | 29 | Intel Skylake | 24 | 192 | High-memory Skylake |
| hpc-n[855-861] | 7 | Intel Cascade Lake | 24 | 128 | AVX-512 |
| hpc-n[895-897] | 3 | Intel Cascade Lake | 16 | 128 | Lower-core nodes |
| hpc-n[894,898-905,916-923] | 17 | Intel Cascade Lake | 32 | 128 | Standard Cascade |
| hpc-n[892-893,876] | 3 | Intel Cascade Lake | 32 | 192 | Higher-memory |
| hpc-n[907-911,913,927-930] | 10 | Intel Cascade Lake | 32 | 370 | Large-memory nodes |
| hpc-n[864-867] | 4 | Intel Cascade Lake | 32 | 384 | Very large-memory |
| hpc-n[883-886] | 4 | AMD EPYC (Rome / epyc2) | 32 | 128 | AVX2 |
| hpc-n[935-942,944-951,953-958,960-966] | 29 | AMD EPYC (Genoa / epyc4) | 32 | 512 | Newest, very high memory |
| hpc-n[975-976] | 2 | AMD EPYC (Genoa / epyc4) | 48 | 128 | Higher core count |
| hpc-n[969-970] | 2 | AMD EPYC (Genoa / epyc4) | 64 | 512 | Highest core count + RAM |
| hpc-n[971-974] | 4 | Intel (64-core AVX-512 class) | 64 | 256 | High-core Intel nodes |

---

$$

53 + 12 + 168 + 29 + 7 + 3 + 17 + 3 + 10 + 4 + 4 + 29 + 2 + 2 + 4 = 347 \text{ nodes}
$$

---

<aside>

### **`grant` Partition — Quick Hardware Snapshot**

- **Total Nodes:** 234
- **Standard Tier (Skylake, 24 cores, ~92 GB):** 146 nodes → default landing spot
- **Medium Tier (Skylake, 24 cores, ~187 GB):** 22 nodes → memory-heavy jobs
- **High Tier (Broadwell, 28 cores, ~251 GB):** 12 nodes → larger RAM fallback
- **Legacy Tier (Broadwell, 28 cores, ~125 GB):** 24 nodes
- **XL Tier (AMD EPYC, 32 cores, 500 GB):** 29 nodes → very large simulations
- **Special:** 1 Cascade Lake node (32 cores, ~185 GB)

### **`public` Partition — Quick Hardware Snapshot**

- **Total Nodes:** 347
- **Standard Tier (Skylake, 24 cores, ~96 GB):** 168 nodes → default landing spot
- **High-Memory Tier (Skylake, 24 cores, ~192 GB):** 29 nodes → memory-heavy CPU jobs
- **Legacy Tier (Broadwell, 28 cores, 128 GB):** 53 nodes → older AVX2 nodes
- **High-Memory Legacy (Broadwell, 28 cores, ~256 GB):** 12 nodes → older but roomy
- **Cascade Lake (mixed):** 44 nodes total
    - 3× (16 cores, 128 GB)
    - 7× (24 cores, 128 GB)
    - 17× (32 cores, 128 GB)
    - 3× (32 cores, ~192 GB)
    - 10× (32 cores, ~370 GB)
    - 4× (32 cores, ~384 GB)
- **EPYC Tier:**
    - **EPYC Genoa/epyc4 (32 cores, 512 GB):** 29 nodes → very large simulations
    - **EPYC Genoa/epyc4 (64 cores, 512 GB):** 2 nodes → max CPU+RAM
    - **EPYC Genoa/epyc4 (48 cores, 128 GB):** 2 nodes → high core count
    - **EPYC Rome/epyc2 (32 cores, 128 GB):** 4 nodes → AVX2
- **Special:** Intel 64-core AVX-512 class (256 GB): 4 nodes

### **🧠 Practical Rule**

- Most jobs → land on **24-core Skylake (~92 GB)**
- Need ~200 GB RAM → target **192 GB Skylake nodes**
- Need >300 GB RAM → must land on **EPYC (500 GB)** nodes

### **🧠 Practical Rule**

- Most jobs → land on **24-core Skylake (~96 GB)**
- Need ~200 GB RAM → target **Skylake 192 GB** nodes
- Need **>300 GB RAM** → target **Cascade Lake 370/384 GB** *or* **EPYC 512 GB** nodes
- Need maximum CPU throughput → target **EPYC 64-core (512 GB)** or **Intel 64-core (256 GB)** nodes
</aside>

## **HPC CPU Generations**

- **SLURM COMMANDS**
    1. Allocate a node: `salloc -p grant --constraint=cascadelake -N 1 --time=01:00:00`
    2. Enter the node: `srun --pty bash`
    3. Inspect CPU: `lscpu` or `lscpu | egrep "CPU\(s\)|Thread|Core|Socket|Model name|CPU max MHz"`
    4. Copy output: `Paste into ChatGPT.`
    5. Get interpreted summary

| **Generation (era)** | **Vendor** | **CPU model (your nodes)** | **Sockets × Cores** | **SMT** | **Key ISA / Features** | **FP64 Peak / Core (base)** | **FP64 Peak / Node (base)** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **Broadwell** (≈2014–2016) | Intel | **Xeon E5-2680 v4 @ 2.40 GHz** | 2 × 14 = 28 | OFF (1) | **AVX2** (no AVX-512), 2 NUMA | 16\times2.4=38.4 GF | 38.4\times28\approx1.08 **TF** |
| **Skylake-SP** (≈2017) | Intel | **Xeon Gold 6126 @ 2.60 GHz** | 2 × 12 = 24 | OFF (1) | **AVX-512** (avx512f/dq/cd/bw/vl), 2 NUMA | 32\times2.6=83.2 GF | 83.2\times24\approx2.00 **TF** |
| **Cascade Lake (R)** (≈2019–2020) | Intel | **Xeon Gold 6226R @ 2.90 GHz** | 2 × 16 = 32 | ON (2) | **AVX-512**, 2 NUMA | 32\times2.9=92.8 GF | 92.8\times32\approx2.97 **TF** |
| **EPYC “epyc4” (Zen4 / Genoa)** (≈2023) | AMD | **EPYC 9124 @ 3.713 GHz** | 2 × 16 = 32 | OFF (1) | **AVX-512**, plus avx512_bf16, avx512_vnni, 2 NUMA | 32\times3.713\approx118.8 GF | 118.8\times32\approx3.80 **TF** |

---

# GPU

## Cluster-wide GPU inventory

| **GPU Model** | **Count** | **Memory per GPU** |
| --- | --- | --- |
| NVIDIA L40S | 8 | 48 GB |
| NVIDIA H100 | 12 | 80 GB |
| NVIDIA A100 | 16 | 40 GB |
| NVIDIA A40 | 10 | 48 GB |
| Quadro RTX 6000 | 29 | 22.7 GB |
| Quadro RTX 5000 | 2 | 16.1 GB |
| V100 SXM2 | 12 | 32.5 GB |
| V100 PCIe | 14 | 32.5 GB |
| P100 | 64 | 16.2 GB |
| GTX 1080 Ti | 100 | 11.1 GB |
| K80 | 20 | 11.4 GB |
| K40 | 10 | 11.4 GB |
| K20 | 20 | 4.7 GB |

## **GPU capability cheat sheet**

| **GPU family** | **Examples** | **Architecture** | **Typical “good for”** |
| --- | --- | --- | --- |
| Pascal | P100 | older | legacy CUDA/HPC, no tensor cores |
| Volta | V100 | older | tensor cores (FP16), solid HPC |
| Turing | RTX 5000/6000 | older | visualization + mixed workloads |
| Ampere | A40/A100 | strong | FP16/TF32, big memory (A100) |
| Ada | L40S | strong | inference/graphics + modern CUDA |
| Hopper | H100/H200 | best | bf16/fp8, fastest DL/HPC |

## GPU Partitions

### `publicgpu` : 139 (in total)

| **Node** | **CPU Generation** | **Cores/node** | **GPU model** | **GPUs/node** | **Notes** |
| --- | --- | --- | --- | --- | --- |
| hpc-n739 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n740 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n741 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n742 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n744 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n745 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n746 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n747 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n748 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n749 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n750 | Intel Skylake | 32 | P100 | 4 | 192G, opa, matlab |
| hpc-n862 | Intel Cascade Lake | 24 | V100 | 4 | 192G, opa, matlab, gputc |
| hpc-n863 | Intel Cascade Lake | 24 | V100 | 4 | 192G, opa, matlab, gputc |
| hpc-n868 | Intel Skylake | 12 | RTX 5000 | 2 | 192G, matlab |
| hpc-n870 | Intel Cascade Lake | 24 | V100 | 4 | 192G, opa, matlab, **NVLink**, gputc |
| hpc-n871 | Intel Cascade Lake | 32 | RTX 6000 | 3 | 192G, **ib**, gputc |
| hpc-n872 | Intel Cascade Lake | 32 | RTX 6000 | 2 | 192G, **ib**, gputc |
| hpc-n873 | Intel Cascade Lake | 32 | RTX 6000 | 3 | 192G, **ib**, gputc |
| hpc-n875 | Intel Cascade Lake | 32 | RTX 6000 | 3 | 192G, **ib**, gputc |
| hpc-n878 | Intel Cascade Lake | 24 | V100 | 4 | 192G, **ib**, **NVLink**, gputc |
| hpc-n881 | Intel Cascade Lake | 32 | RTX 6000 | 3 | 192G, **ib**, gputc |
| hpc-n882 | Intel Cascade Lake | 32 | RTX 6000 | 3 | 192G, **ib**, gputc |
| hpc-n888 | AMD EPYC2 (Rome) | 48 | A100 | 2 | 128G, **ib**, **40GB** |
| hpc-n889 | AMD EPYC2 (Rome) | 48 | A100 | 2 | 128G, **ib**, **40GB** |
| hpc-n890 | AMD EPYC2 (Rome) | 48 | A100 | 2 | 128G, **ib**, **40GB** |
| hpc-n891 | AMD EPYC2 (Rome) | 48 | A100 | 2 | 128G, **ib**, **40GB** |
| hpc-n925 | Intel Cascade Lake | 32 | A40 | 2 | 192G |
| hpc-n926 | Intel Cascade Lake | 32 | A40 | 2 | 192G |
| hpc-n931 | AMD EPYC3 (Milan) | 64 | A100 | 4 | 128G, **80GB**, HGX |
| hpc-n932 | AMD EPYC3 (Milan) | 64 | H100 | 4 | **1024G**, **40GB** |
| hpc-n933 | Intel (AVX-512 class) | 32 | H100 | 4 | **512G**, **NVLink**, **80GB** |
| hpc-n934 | Intel (AVX-512 class) | 32 | H100 | 4 | **512G**, **NVLink**, **80GB** |
| hpc-n967 | AMD EPYC3 (Milan) | 64 | L40S | 4 | **512G**, **45GB** |
| hpc-n968 | AMD EPYC3 (Milan) | 64 | L40S | 4 | **512G**, **45GB** |
| hpc-n978 | AMD EPYC4 (Genoa) | 64 | H200 | 8 | **2300G** |
| hpc-n979 | AMD EPYC4 (Genoa) | 64 | H200 | 8 | **2300G** |
| hpc-n981 | AMD EPYC4 (Genoa) | 64 | H200 | 8 | **2300G** |

### `grantgpu` : 50 (in total)

| **Node** | **CPU Generation** | **Cores/node** | **GPU model** | **GPUs/node** | **Notes** |
| --- | --- | --- | --- | --- | --- |
| hpc-n862 | Intel Cascade Lake | 24 | V100 | 4 | 192G, opa, matlab, gputc |
| hpc-n863 | Intel Cascade Lake | 24 | V100 | 4 | 192G, opa, matlab, gputc |
| hpc-n878 | Intel Cascade Lake | 24 | V100 | 4 | 192G, ib, **NVLink**, gputc |
| hpc-n881 | Intel Cascade Lake | 32 | RTX 6000 | 3 | 192G, ib, gputc |
| hpc-n882 | Intel Cascade Lake | 32 | RTX 6000 | 3 | 192G, ib, gputc |
| hpc-n888 | AMD EPYC2 (Rome) | 48 | A100 | 2 | 128G, ib, **40GB** |
| hpc-n889 | AMD EPYC2 (Rome) | 48 | A100 | 2 | 128G, ib, **40GB** |
| hpc-n890 | AMD EPYC2 (Rome) | 48 | A100 | 2 | 128G, ib, **40GB** |
| hpc-n891 | AMD EPYC2 (Rome) | 48 | A100 | 2 | 128G, ib, **40GB** |
| hpc-n978 | AMD EPYC4 (Genoa) | 64 | H200 | 8 | **2300G** |
| hpc-n979 | AMD EPYC4 (Genoa) | 64 | H200 | 8 | **2300G** |
| hpc-n981 | AMD EPYC4 (Genoa) | 64 | H200 | 8 | **2300** |