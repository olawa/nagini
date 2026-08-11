__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"

"""Generate Hydra-style copy rules for output entries with an input path."""


output_directory = pathlib.Path(output_spec.get("directory", "./"))
copy_rule_strings = []
for file_definition in output_spec["files"]:
    if file_definition.get("when") and not config.get(file_definition["when"], False):
        continue
    if file_definition["input"] is None:
        continue

    rule_name = "_copy_{}".format(
        "_".join(re.split(r"\W{1,}", file_definition["name"].strip().lower()))
    )
    input_file = pathlib.Path(file_definition["input"])
    output_file = output_directory / pathlib.Path(file_definition["output"])
    copy_resources = config.get("_copy", {})
    mem_mb = copy_resources.get("mem_mb", config["default_resources"]["mem_mb"])
    mem_per_cpu = copy_resources.get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
    partition = copy_resources.get("partition", config["default_resources"]["partition"])
    threads = copy_resources.get("threads", config["default_resources"]["threads"])
    time = copy_resources.get("time", config["default_resources"]["time"])
    copy_container = copy_resources.get("container", config["default_container"])

    copy_rule_strings.append(
        "\n".join(
            [
                f'@workflow.rule(name="{rule_name}")',
                f'@workflow.input("{input_file}")',
                f'@workflow.output("{output_file}")',
                f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
                f'@workflow.container("{copy_container}")',
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb="{mem_mb}", '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
                f'@workflow.shellcmd("{copy_container}")',
                "@workflow.run\n",
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, "
                "log, version, rule, conda_env, container_img, singularity_args, use_singularity, "
                "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
                "shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):",
                '\tshell("(cp --preserve=timestamps {input[0]} {output[0]}) &> {log}", '
                "bench_record=bench_record, bench_iteration=bench_iteration)\n\n",
            ]
        )
    )

if copy_rule_strings:
    exec(compile("\n".join(copy_rule_strings), "copy_result_files", "exec"), workflow.globals)
