---
name: riboclaw-visualization
description: Produce RiboClaw-friendly structured process visualization output (chain + step cards + image artifacts). Use when the user asks to visualize reasoning, show processing DAG/chain, render step-by-step execution, or drive the RiboClaw frontend timeline.
---

# riboclaw-visualization

When running tasks for the RiboClaw UI, include a final JSON code block for UI rendering.

## Required JSON Schema

```json
{
  "title": "任务标题",
  "summary": "一句话总结",
  "chain": [
    {"id":"NODE_01","module":"模块名","status":"success|running|failed","meta":"简短元数据"}
  ],
  "steps": [
    {
      "node_id":"NODE_01",
      "module":"模块名",
      "report":"用自然语言自由汇报这一步做了什么、发现了什么、关键数值和结论。像写实验笔记一样，直接说重点。",
      "images":["图片产物的绝对路径，可为空数组"]
    }
  ]
}
```

## Rules

1. Ensure chain order matches real execution order.
2. Include `failed` status when any critical step fails.
3. Put generated plot/image paths into `images`.
4. Do not fabricate files; only include paths that actually exist.
5. `report` field should be natural, conversational — describe findings and key numbers directly, not a rigid template.
