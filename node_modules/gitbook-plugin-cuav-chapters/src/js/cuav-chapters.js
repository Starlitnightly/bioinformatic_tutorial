require(['gitbook', 'jQuery'], function(gitbook, $) {

  var TOGGLE_CLASSNAME = 'cuav-expanded';
  var CHAPTER_CLASS = '.chapter';
  var CHAPTER_CLASSNAME = 'chapter';
  var ARTICLES_CLASS = '.articles';
  var ARTICLES_CLASSNAME = 'articles';
  var TRIGGER_TEMPLATE = '<i class="cuav-trigger fa"></i>';

  // 展开目录的最大个数，当目录数小于此数值时，自动张开全部目录
  var summaryMaxSize = 20;

  // 统计目录数
  var chapterSize = 0;

  /**
   * 初始化
   * @function
   */
  var init = function() {
    chapterSize = 0;

    var cuavChaptersConfig
      = gitbook.state.config.pluginsConfig['cuav-chapters'];
    
    if (cuavChaptersConfig.summaryMaxSize) { // 获取设置的目录展开最大值
      summaryMaxSize = cuavChaptersConfig.summaryMaxSize;
    }

    if (cuavChaptersConfig.useLimitExpanded) { // 使用限制展开、收起目录
      limitExpanded();
    } else {
      bindEvent();
      expand($('li.chapter.active'));
    }

    ajaxSummary(chapterSize);

  }

  /**
   * 绑定事件
   * @function
   */
  var bindEvent = function() {
    var $chapterLinkNode
      = $(ARTICLES_CLASS).parent(CHAPTER_CLASS).children('a,span');

    for (var i = 0; i < $chapterLinkNode.length; ++i) {
      bindClickEvent($chapterLinkNode.eq(i), $(TRIGGER_TEMPLATE));
    }
  }

  /**
   * 绑定目录点击事件
   * @function
   * @param  {jQuery} $node  目录节点（<a>/<span>）
   * @param  {jQuery} $iNode 目录的小图标
   */
  var bindClickEvent = function($node, $iNode) {
    
    if ($node.is('span')
      || ($node.is('a') 
        && $node.attr('href') 
        && ($node.attr('href') === 'javascript:;'
          || $node.attr('href') === '#'))) { // 目录点击展开、收起
      $node.on('click', function(event) {
        event.stopPropagation();
        toggle($(event.target).parent(CHAPTER_CLASS));
      });
    }

    $iNode.on('click', function(event) { // 处理点击目录图标
      event.preventDefault();
      event.stopPropagation();
      toggle($(event.target).parent().parent(CHAPTER_CLASS));
    });

    $node.append($iNode);
  }

  /**
   * 展开、收起目录
   * @function
   * @param  {jQuery} $chapter 目录节点
   */
  var toggle = function($chapter) {
    if ($chapter.hasClass(TOGGLE_CLASSNAME)) { // 收起目录
      $chapter.removeClass(TOGGLE_CLASSNAME);
    } else { // 展开目录
      $chapter.addClass(TOGGLE_CLASSNAME);
    }
  }

  /**
   * 展开当前目录及其父目录
   * @function
   * @param  {jQuery} $chapter 要展开的目录
   */
  var expand = function($chapter) {
    $chapter.addClass(TOGGLE_CLASSNAME);

    var $chapterParent = $chapter.parent();
    if (!$chapterParent.hasClass('summary')
      && !$chapterParent.hasClass('book-summary')
      && $chapter.length != 0) { // 递归展开父目录
      expand($chapterParent);
    }
  }

  /**
   * 使用限制展开目录
   * @function
   */
  var limitExpanded = function() {

    var calcSize = function(item) { // 计算目录数
      var $item = $(item);
      if ($item.is('ul')) {
        $chapters = $item.children('li.chapter');
        chapterSize += $chapters.length;
        $.each($chapters, function(index, val) {
           calcSize(val);
        });

      } else if ($item.is('li')) {
        if ($item.has('ul.articles')) {
          calcSize($item.children('ul.articles'));
        }
      }

    }

    calcSize('.summary');

    if (chapterSize >= summaryMaxSize) { // 收起目录
      bindEvent();
      expand($('li.chapter.active')); 
    } else { // 展开目录
      $.each($('.chapter > .articles'), function(index, val) {
        $(val).css('max-height', '9999px');
      });
    }
  }

  /**
   * 远程目录
   * @function
   */
  var ajaxSummary = function(chapterSize) {

    var cuavChaptersConfig
      = gitbook.state.config.pluginsConfig['cuav-chapters'];

    if (cuavChaptersConfig.navUrl) { // 有配置远程配置目录 url，添加远程配置目录
      $.ajax(cuavChaptersConfig.navUrl, {
        type: 'get',
        dataType: 'json',
        cache: false,
        success: function(data) {
          if (data) {
            var newChapterSize = chapterSize;

            for (var i = 0; i < data.length; ++i) {
              ++newChapterSize;
              var item = data[i];
              if (item.links) {
                newChapterSize += item.links.length;
              }
            }

            buildSummary(data, newChapterSize);
          }
        }
      });
    }

    /**
     * 构建远程配置目录
     * @function
     * @param  {array} items            远程目录
     * @param  {number} newChapterSize 旧目录个数
     */
    var buildSummary = function(items, newChapterSize) { // 构建目录
      var $summary = $('.summary');

      // 决定构建后的整体目录是何种状态：
      // -1：之前展开，构建后展开
      // 0：之前展开，构建后收起
      // 1：之前收起，构建后收起
      var status;

      var cuavChaptersConfig
        = gitbook.state.config.pluginsConfig['cuav-chapters'];

      if (cuavChaptersConfig.useLimitExpanded) {
        if (chapterSize >= summaryMaxSize) {
          status = 1;
        } else {
          if (newChapterSize >= summaryMaxSize) {
            status = 0;
          } else {
            status = -1;
          }
        }
      } else {
        status = 1;
      }

      for (var i = 0; i < items.length; ++i) { // 循环创建目录
        var item = items[i];
        if (item.name) { // 如果没有目录名，则不构建此目录
          var $chapterNode = $('<li>').addClass(CHAPTER_CLASSNAME);
          var $chapterLinkNode
            = $('<a>')
                .attr('href', item.url ? item.url : 'javascript:;')
                .append(item.name);
          $chapterNode.append($chapterLinkNode);


          if (item.links) { // 有二级目录
            if (status === 1) { // 绑定目录点击事件
              bindClickEvent($chapterLinkNode, $(TRIGGER_TEMPLATE));
            }

            var links = item.links;
            var $articlesNode = $('<ul>').addClass(ARTICLES_CLASSNAME);

            if (status === -1) {
              $articlesNode.css("max-height", "9999px");
            }

            for (var j = 0; j < links.length; ++j) { // 循环构建二级目录
              var link = links[j];
              if (link.name) {
                var $chapterSecondNode = $('<li>').addClass(CHAPTER_CLASSNAME);
                var $chapterSecondLinkNode
                  = $('<a>')
                    .attr('href', link.url ? link.url : 'javascript:;')
                    .append(link.name);

                $articlesNode.append(
                  $chapterSecondNode.append($chapterSecondLinkNode));
              }

              $chapterNode.append($articlesNode);
            }
          }

          $summary.append($chapterNode);
        }
      }

      if (status === 0) { // 之前展开，构建后收起
        $.each($('.chapter > .articles'), function(index, val) {
          $(val).css('max-height', '');
        });
        bindEvent();
        expand($('li.chapter.active'));
      }

    }
  }

  gitbook.events.bind('page.change', function() {
    init();
  });
});