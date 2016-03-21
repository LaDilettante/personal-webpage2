---
layout: page
title: "Statistics Book"
description: "My notes when 'Exercise is left to the reader'"
header-img: "img/about-bg.jpg"
---

<ul>
{% for chap in site.statbook %}
  <li>
    {% if chap.title %}<a href="{{ chap.url | prepend: site.baseurl }}">{{ chap.title }}</a>{% endif %}
  </li>
{% endfor %}
</ul>